#pragma once
#include <iostream>
#include <utility>
#include <boost/any.hpp>
#include <boost/unordered_map.hpp>
#include <vector>
#include <map>
#include <string>
#include <functional>
#include <typeinfo>
#include <typeindex>
#include <boost/optional.hpp>
#ifndef __CUDACC__
#include <filesystem>
#include <variant>
#endif
//#include <gvdb.h>
#include <eigen3/Eigen/Dense>

// #if !__has_include(<cuda_runtime.h>)
using float2 = Eigen::Vector2f;
using float3 = Eigen::Vector3f;
using float4 = Eigen::Vector4f;
using  int2 =Eigen::Vector2i;
using  int3 =Eigen::Vector3i;
using  int4 =Eigen::Vector4i;
using double2 =Eigen::Vector2d;
using double3 =Eigen::Vector3d;
using double4 =Eigen::Vector4d;
// #endif

namespace YAML {
	class Node;
}

namespace detail {
	template<typename Test, template<typename...> class Ref>
	struct is_specialization : std::false_type {};
	template<template<typename...> class Ref, typename... Args>
	struct is_specialization<Ref<Args...>, Ref> : std::true_type {};
	struct iAny : boost::any {
		using boost::any::any;
		template<class T> operator T() const { return boost::any_cast<T>(*this); }
		template<class T> operator T& () { return boost::any_cast<T&>(*this); }
	};
	struct varAny {
		boost::optional<iAny> val;
		boost::optional<iAny> valVec;
		bool isVec = false;

		template<typename T>
		varAny(T var) :val(var), isVec(false) {}
		template<typename T>
		varAny(std::vector<T> var) : valVec(var), isVec(true) {}

		template<class T> operator std::vector<T>() const {
			return boost::any_cast<std::vector<T>>(valVec.value());
		}
		template<class T> operator std::vector<T>& () {
			return boost::any_cast<std::vector<T>&>(valVec.value());
		}
		template<class T> operator T() const { return boost::any_cast<T>(val.value()); }
		template<class T> operator T& () { return boost::any_cast<T&>(val.value()); }
	};
}
struct Range { detail::iAny min, max;};
struct ParameterSettings {
	std::string description = "";
	bool constant = false;
	bool hidden = false;
	std::vector<std::type_index> alternativeTypes{};
	boost::optional<Range> range;
	std::vector<detail::iAny> presets{};
};
struct Parameter {
	std::string identifier;
	std::string identifierNamespace;

	detail::varAny param;
	std::type_index type;
	ParameterSettings properties;
};


std::pair<std::string, std::string> split(std::string s);
class ParameterManager {

	bool parameterExists(std::string s);
	std::string resolveParameter(std::string s);
	void parseTree(YAML::Node root);
public:
	bool isAmbiguous(std::string ident);

	std::map<std::string, Parameter*> parameterList;
	std::multimap<std::string, std::string> qualifiedIdentifierList;
	boost::unordered_map<std::type_index, std::function<detail::iAny(const YAML::Node&)>> decoders;
	boost::unordered_map<std::type_index, std::function<YAML::Node(const detail::iAny&)>> encoders;
	boost::unordered_map<std::type_index, std::function<void(Parameter&)>> uiFunctions;

	ParameterManager();
	void addDecoder(std::type_index ty, std::function<detail::iAny(const YAML::Node&)> fn);
	void addEncoder(std::type_index ty, std::function<YAML::Node(const detail::iAny&)> fn);
	void addUifunction(std::type_index ty, std::function<void(Parameter&)> fn);
	// static ParameterManager& instance();
	template<typename T>
	void newParameter(std::string identifier, T defaultValue) {
		if (parameterList.find(identifier) != parameterList.end())
			throw std::invalid_argument("Parameter " + identifier + " already exists.");
		auto nsid = split(identifier);
		auto ns = nsid.first;
		auto id = nsid.second;
		qualifiedIdentifierList.insert(std::make_pair(id, identifier));
		parameterList[identifier] = new Parameter{ id, ns, defaultValue, typeid(T) };
	}
	template<typename T>
	void newParameter(std::string identifier, T defaultValue, ParameterSettings&& props) {
		if (parameterList.find(identifier) != parameterList.end())
			throw std::invalid_argument("Parameter " + identifier + " already exists.");
		auto nsid = split(identifier);
		auto ns = nsid.first;
		auto id = nsid.second;
		qualifiedIdentifierList.insert(std::make_pair(id, identifier));
		parameterList[identifier] = new Parameter{ id, ns, defaultValue, typeid(T), std::move(props) };
	}
	Parameter& getParameter(std::string identifier);
	detail::varAny& get(std::string identifier);
	template<typename T>
	std::enable_if_t<!detail::is_specialization<T, std::vector>::value, T&> get(std::string identifier) {
		return boost::any_cast<T&>(getParameter(identifier).param.val.value());
	}
	template<typename T>
	std::enable_if_t<detail::is_specialization<T, std::vector>::value, std::vector<typename T::value_type>&> get(std::string identifier) {
		return boost::any_cast<std::vector<typename T::value_type>&>(getParameter(identifier).param.valVec.value());
	}



	void load(std::string filename);
	void loadDirect(std::string yaml);
	std::string buildTree();
	void loadTree(YAML::Node);
#ifndef NO_GUI
	void buildImguiWindow(bool* p_open);
#endif
	void init();
};


#define callEncoder(ty, x) pm.encoders[typeid(ty)](x);
#define callDecoder(ty, n, alt) n ? boost::any_cast<ty>(pm.decoders[typeid(ty)](n)) : alt;
#define callUI(ty, variable, member)\
{\
	auto xParam = Parameter{ parameter.identifier + member, parameter.identifierNamespace, variable, typeid(ty), parameter.properties };\
	pm.uiFunctions[typeid(ty)](xParam);\
	if(!parameter.properties.constant)\
		variable = boost::any_cast<ty>(xParam.param.val.value());\
}
#define customVector(ty)\
	pm.addEncoder(typeid(std::vector<ty>), [this](const detail::iAny& any) {\
	printf(".");\
		auto var = boost::any_cast<std::vector<ty>>(any);\
		printf("-");\
		auto node = YAML::Node();\
		printf(",");\
		for (const auto& elem : var){\
			printf("=");\
			detail::iAny any = detail::iAny(elem);\
			node.push_back(pm.encoders[typeid(std::decay_t<ty>)](any));\
		}\
		printf(":");\
		return node;\
		});\
	pm.addDecoder(typeid(std::vector<ty>), [this](const YAML::Node& node) {\
		std::vector<ty> vec;\
		for (const auto& child : node)\
			vec.push_back((ty)pm.decoders[typeid(ty)](child));\
		return detail::iAny(vec);\
		}); //\
	// pm.addUifunction(typeid(std::vector<ty>), [this](Parameter& parameter) {\
	// 	std::vector<ty>& vec = boost::any_cast<std::vector<ty>&>(parameter.param.valVec.value());\
	// 	if (parameter.properties.hidden) return;\
	// 	ImGui::PushID(parameter.identifier.c_str());\
	// 	ImGui::Text(parameter.identifier.c_str());\
	// 	if (!parameter.properties.constant) {\
	// 		ImGui::SameLine();\
	// 		if (ImGui::Button("+"))\
	// 			vec.push_back(ty());\
	// 		ImGui::SameLine();\
	// 		if (ImGui::Button("-"))\
	// 			vec.pop_back();\
	// 	}\
	// 	int32_t i = 0;\
	// 	for (auto& elem : vec) {\
	// 		ImGui::Indent();\
	// 		ImGui::PushID(i);\
	// 		auto eParam = Parameter{ parameter.identifier + "[" + std::to_string(i++) + "]", parameter.identifierNamespace, elem, typeid(ty), parameter.properties };\
	// 		pm.uiFunctions[typeid(ty)](eParam);\
	// 		ImGui::Unindent();\
	// 		ImGui::PopID();\
	// 		if (!parameter.properties.constant)\
	// 			elem = (ty) eParam.param.val.value();\
	// 	}\
	// 	ImGui::PopID();\
	// 	});
