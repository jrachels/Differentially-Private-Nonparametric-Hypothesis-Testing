#include "DifPrivNonparametricHypthTestingTests/catch.hpp"

#include "DifferentiallyPrivateNonparametricHypthothesisTesting/Include/dpnht.h"

TEST_CASE("Basic run_time test for Private Kruskal Wallis") {

	std::random_device  rd{};
	std::mt19937 mt(rd());

	std::vector<std::vector<double>> raw_data;
	for (int i = 0; i < 7; i++) {
		std::vector<double> column;

		for (int j = 0; j < i; j++) {
			column.push_back(i * j);
		}

		raw_data.push_back(column);
	}

	auto result1 = dpnht::PrivateKruskalWallis<double, std::mt19937>(raw_data, .1, mt, 20, false);
	auto result2 = dpnht::PrivateKruskalWallis<double>(raw_data, .1, 20, false);

	auto result3 = dpnht::PrivateKruskalWallis<double, std::mt19937>(raw_data, .1, mt, 20, true);
	auto result4 = dpnht::PrivateKruskalWallis<double>(raw_data, .1, 20, true);

}