#include "DifPrivNonparametricHypthTestingTests/catch.hpp"

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/Include/dpnht.h"

TEST_CASE("Basic run_time test for Private Wilcoxon") {

	std::random_device  rd{};
	std::mt19937 mt(rd());

	std::vector<std::pair<double, double>> raw_data3(50);

	for (int i = 0; i < 50; i++) {
		raw_data3[i].first = (7 * i - 50);
		raw_data3[i].first = (4 * i - 12);
	}

	auto result7 = dpnht::PrivateWilcoxon(raw_data3, .1, mt, 24);
	auto result8 = dpnht::PrivateWilcoxon(raw_data3, .1, 24);

}