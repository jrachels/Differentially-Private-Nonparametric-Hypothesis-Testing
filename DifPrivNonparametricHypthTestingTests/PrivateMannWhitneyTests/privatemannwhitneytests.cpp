#include "DifPrivNonparametricHypthTestingTests/catch.hpp"

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/Include/dpnht.h"

TEST_CASE("Basic run_time test for Private Mann Whitney") {
	std::random_device  rd{};
	std::mt19937 mt(rd());


	std::pair<std::vector<double>, std::vector<double>> raw_data2;

	for (int j = 0; j < 12; j++) {
		raw_data2.first.push_back(5 * j);
	}
	for (int j = 0; j < 17; j++) {
		raw_data2.second.push_back(4 * j);
	}


	auto result5 = dpnht::PrivateMannWhitney(raw_data2, .1, .3, .2, mt, 24);
	auto result6 = dpnht::PrivateMannWhitney(raw_data2, .1, .3, .2, 24);

}