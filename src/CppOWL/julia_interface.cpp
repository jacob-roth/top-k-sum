#include <jlcxx/jlcxx.hpp>
#include "project_to_OWL_ball.h"  // Header that declares evaluateProx and other functions

// Define the module's interface for Julia
JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.method("evaluateProx", [](std::vector<double>& z_in, std::vector<double>& w, double epsilon, std::vector<double>& x_out, bool sorted_and_positive) {
        evaluateProx(z_in, w, epsilon, x_out, sorted_and_positive);
    });
}
