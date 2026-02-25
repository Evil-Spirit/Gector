#pragma once
#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Minimal test framework
// ---------------------------------------------------------------------------
struct Test {
    std::string              name;
    std::function<void()>    fn;
};

inline std::vector<Test>& tests() {
    static std::vector<Test> g;
    return g;
}

inline void registerTest(std::string name, std::function<void()> fn) {
    tests().push_back({std::move(name), std::move(fn)});
}

// ---------------------------------------------------------------------------
// Assertion macros
// ---------------------------------------------------------------------------
#define ASSERT_TRUE(expr) \
    do { if (!(expr)) throw std::runtime_error("ASSERT_TRUE failed: " #expr); } while(0)

#define ASSERT_FALSE(expr) \
    do { if ((expr))  throw std::runtime_error("ASSERT_FALSE failed: " #expr); } while(0)

#define ASSERT_EQ(a, b) \
    do { if (!((a) == (b))) throw std::runtime_error("ASSERT_EQ failed: " #a " != " #b); } while(0)

#define ASSERT_NEAR(a, b, eps) \
    do { \
        double _a = static_cast<double>(a), _b = static_cast<double>(b); \
        if (std::abs(_a - _b) >= static_cast<double>(eps)) { \
            throw std::runtime_error( \
                std::string("ASSERT_NEAR failed: |") + std::to_string(_a) + \
                " - " + std::to_string(_b) + "| >= " + std::to_string(eps)); \
        } \
    } while(0)

#define ASSERT_THROWS(expr) \
    do { \
        bool _threw = false; \
        try { (void)(expr); } catch (...) { _threw = true; } \
        if (!_threw) throw std::runtime_error("ASSERT_THROWS: no exception thrown for " #expr); \
    } while(0)

// ---------------------------------------------------------------------------
// Registration helper (used at file scope)
// ---------------------------------------------------------------------------
struct TestRegistrar {
    TestRegistrar(const char* name, std::function<void()> fn) {
        registerTest(name, std::move(fn));
    }
};

#define TEST(suite, name) \
    static void test_##suite##_##name(); \
    static TestRegistrar reg_##suite##_##name(#suite "/" #name, test_##suite##_##name); \
    static void test_##suite##_##name()
