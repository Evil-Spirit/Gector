#include "test_framework.h"
#include <iostream>

int main() {
    int passed = 0, failed = 0;
    for (const auto& t : tests()) {
        try {
            t.fn();
            std::cout << "[PASS] " << t.name << "\n";
            ++passed;
        } catch (const std::exception& e) {
            std::cout << "[FAIL] " << t.name << ": " << e.what() << "\n";
            ++failed;
        } catch (...) {
            std::cout << "[FAIL] " << t.name << ": unknown exception\n";
            ++failed;
        }
    }
    std::cout << "\n" << passed + failed << " tests: "
              << passed << " passed, " << failed << " failed\n";
    return failed > 0 ? 1 : 0;
}
