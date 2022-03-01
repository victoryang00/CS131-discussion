#include <iostream>
//#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <assert.h>
#include <utility>
#include <benchmark/benchmark.h>

using namespace std;
//using namespace cv;

random_device seed;
mt19937 gen(seed());
const int x = 19;
const int y = 9;
constexpr const char str1[] = "yangyw";
constexpr const char str2[] = "ISvegetable";

// Intro to lambda programming
//static void test_lambda(benchmark::State &state) {
//    for (auto _: state) {
//        vector<Point2i> vec;
//        vec.push_back(Point2i(3, 1));
//        vec.push_back(Point2i(3, 3));
//        vec.push_back(Point2i(2, 3));
//        vec.push_back(Point2i(2, 1));
//        vec.push_back(Point2i(1, 3));
//        vec.push_back(Point2i(1, 1));
//        vec.push_back(Point2i(2, 2));
//        vec.push_back(Point2i(1, 2));
//
//        sort(vec.begin(), vec.end(),
//             [](auto pt1, auto pt2) { return (pt1.x != pt2.x) ? (pt1.x < pt2.x) : (pt1.y < pt2.y); });
//    }
//}

//helper function
//bool cmp(Point2i pt1, Point2i pt2) {
//    return (pt1.x != pt2.x) ? (pt1.x < pt2.x) : (pt1.y < pt2.y);
//}

//static void test_without_lambda(benchmark::State &state) {
//    for (auto _: state) {
//        vector<Point2i> vec;
//        vec.push_back(Point2i(3, 1));
//        vec.push_back(Point2i(3, 3));
//        vec.push_back(Point2i(2, 3));
//        vec.push_back(Point2i(2, 1));
//        vec.push_back(Point2i(1, 3));
//        vec.push_back(Point2i(1, 1));
//        vec.push_back(Point2i(2, 2));
//        vec.push_back(Point2i(1, 2));
//
//        sort(vec.begin(), vec.end(), cmp);
//    }
//}

// Intro to copy Iterator
static void test_copy_iterator(benchmark::State &state) {
    for (auto _:state) {
        vector<int> nums;
        nums.push_back(3);
        nums.push_back(9);
        nums.push_back(55);
        nums.push_back(4);
        //output:3,9,55,4
        sort(nums.begin(), nums.end());
//        copy(begin(nums), end(nums) );
        //output:3,4,9,55
        sort(nums.begin(), nums.end(), greater<>());
//        copy(begin(nums), end(nums)、 );
        //output:55,9,4,3
    }
};

//  mul 3 5 -> 3 * 5
template<int a, int b>
struct mul {
    static constexpr int value = a * b;
};

static void test_mul(benchmark::State &state) {
    while (state.KeepRunning()) {
        auto r = mul<x, y>::value;
        assert(r == x * y);
    }
}


// split int-> ‘int’
template<const char *str, int len, char... suffix>
struct append {
    static constexpr const char *value() {
        return append<str, len - 1, str[len - 1], suffix...>::value();
    }
};

template<const char *str, char... suffix>
struct append<str, 0, suffix...> {
    static const char value_str[];

    static constexpr const char *value() {
        return value_str;
    }
};

template<const char *str, char... suffix>
const char append<str, 0, suffix...>::value_str[] = {suffix..., 0};

template<typename T>
struct base_typename_struct;

template<>
struct base_typename_struct<int> {
    static constexpr const char name[] = "int";
};

template<typename T, char... suffix>
struct typename_struct {
    typedef base_typename_struct<T> base;

    static const char *name() {
        return append<base::name, sizeof(base::name) - 1, suffix...>::value();
    }
};

template<typename T, char... suffix>
struct typename_struct<T *, suffix...> :
        public typename_struct<T, '*', suffix...> {
};

static void test_split_type(benchmark::State &state) {
    while (state.KeepRunning()) {
        auto r = typename_struct<int ****>::name();
        assert(r[0] == 'i');
    }
}

// con ‘yangyw’‘ISvegetable’-> ‘yangywISvegetable’
template<int...I> using is = std::integer_sequence<int, I...>;
template<int N> using make_is = std::make_integer_sequence<int, N>;

constexpr auto size(const char *s) {
    int i = 0;
    while (*s != 0) {
        ++i;
        ++s;
    }
    return i;
}

template<const char *, typename, const char *, typename>
struct concat_impl;

template<const char *S1, int... I1, const char *S2, int... I2>
struct concat_impl<S1, is<I1...>, S2, is<I2...>> {
    static constexpr const char value[]{S1[I1]..., S2[I2]..., 0};
};

template<const char *S1, const char *S2>
constexpr auto concat{
        concat_impl<S1, make_is<size(S1)>, S2, make_is<size(S2)>>::value
};

static void test_con(benchmark::State &state) {
    while (state.KeepRunning()) {
        auto r = concat<str1, str2>;
        assert(r[0] == 'y');
    }
}

// Factorial
template<long long N>
struct Factorial {
    enum {
        value = N * Factorial<N - 1>::value
    };
};

template<>
struct Factorial<0> {
    enum {
        value = 1
    };
};

// pattern matching
struct Zero {
    enum {
        value = 0
    };
};

template<typename N>
struct Succ {
    enum {
        value = N::value + 1
    };
};

template<typename N>
struct MatchOne {
    enum {
        value = 0
    };
};

template<>
struct MatchOne<Succ<Zero> > {
    enum {
        value = 1
    };
};

static void test_fatorial(benchmark::State &state) {
    for (auto _: state) {
        Factorial<20>::value;
    }
}

static void test_match(benchmark::State &state) {
    for (auto _:state) {
        MatchOne<Zero>::value;

        MatchOne<Succ<Zero> >::value;

        MatchOne<Succ<Succ<Zero> > >::value;
    }
}

// palindrome

/** palindrome [] = True
palindrome [_] = True
palindrome (elem:rest) = (elem == last rest) && (palindrome(init rest))
*/
template<class T>
bool palindrome(T &v) {
    return std::equal(v.begin(), v.end(), v.rbegin());
}

static void test_palindrome(benchmark::State &state) {
    for (auto _:state) {
        int ar[] = {1, 2, 3, 2, 1};
        vector<int> v(ar, ar + 5);
        assert(palindrome(v) == true);
    }
}

//branching 
template<bool Condition, typename TrueResult, typename FalseResult>
class if_;

template<typename TrueResult, typename FalseResult>
struct if_<true, TrueResult, FalseResult> {
    typedef TrueResult result;
};

template<typename TrueResult, typename FalseResult>
struct if_<false, TrueResult, FalseResult> {
    typedef FalseResult result;
};

static void test_branching(benchmark::State &state) {
    for (auto _:state) {
        typename if_<true, int, void *>::result number(3);
        typename if_<false, int, void *>::result pointer(&number);

        typedef typename if_<(sizeof(void *) > sizeof(uint32_t)), uint64_t, uint32_t>::result
                integral_ptr_t;

        integral_ptr_t converted_pointer = reinterpret_cast<integral_ptr_t>(pointer);
    }
}

//BENCHMARK(test_lambda);
//BENCHMARK(test_without_lambda);
BENCHMARK(test_fatorial);
BENCHMARK(test_match);
//BENCHMARK(test_copy_iterator);
BENCHMARK(test_mul);
BENCHMARK(test_split_type);
BENCHMARK(test_con);
BENCHMARK(test_palindrome);
BENCHMARK(test_branching);

BENCHMARK_MAIN();
