#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
#include "metric_search.hpp"
// Make fuzzy testing

template<typename T>
struct distance {
    int operator()( const T &lhs,  const T &rhs) const {
	return std::abs(lhs - rhs);
    }
};
struct random_uniform_int {
  //  std::random_device rd;
  std::mt19937 gen;
  std::uniform_int_distribution<int> d;
  random_uniform_int(int min, int max)
      :gen(),d(min,max) {}
  int operator()()  {
    return d(gen);
  }
};
struct random_uniform_real {
  //  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float> d;
  random_uniform_real(float min, float max)
      :gen(),d(min,max) {}
  float operator()()  {
    return d(gen);
  }
};
bool test(int array_size) {
    //    random_uniform_int dgen(std::numeric_limits<float>::min()/2,std::numeric_limits<float>::max()/2);
    random_uniform_real  dgen(-1000000,1000000);
    std::vector<float> data;
    data.reserve(array_size);
    for(int i = 0; i < array_size; i++) {
        data.push_back(dgen());
    }
    metric_search::Tree<float,distance<float>> tr;
    for(auto i : data) {
        tr.insert(i);
        if(!tr.check_covering()) {
            return false;
        }
    }
    for(std::size_t i = 0; i < data.size(); i++) {
        auto root = tr.get_root();
        bool b = root->children.empty();
        tr.erase(root->data);
        if(!b) {
            if(!tr.check_covering()) {
                return false;
            }
        }
    }
    return true;
}

int main(int argc, char ** argv) {
    
    int iterations = 1;
    if(argc == 2)
        iterations = std::stoi(argv[1]);
    random_uniform_int len_gen(1,1000);
    for(int i = 0; i < iterations; i++) {
        int array_size = len_gen();
        bool result = test(array_size);
        if(!result) {
            std::cout << "Error!!! -- " << i << "," << array_size << std::endl;
        }
        std::cout << "Iteration #" << i << std::endl;
    }
    return 0;
}
