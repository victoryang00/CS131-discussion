//safe and sound

#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/foreach.hpp>


#include <iomanip>
#include <cmath>

double binom(int n, int k) { return 1 / ((n + 1) * std::beta(n - k + 1, k + 1)); }

class STD_Suggestion {
public:
    STD_Suggestion() : STD_Suggestion{false, false, false} {}

    STD_Suggestion(bool safe, bool modern, bool free)
            : SMF{
            {"safe", safe},
            {"modern", modern},
            {"free", free}
    } {
        suggestion = "suggest";
    }

    void print_suggestion();

private:
    void set_suggestion();

    std::unordered_map<std::string, bool> SMF;
    std::string suggestion;
};

void STD_Suggestion::print_suggestion() {
    set_suggestion();
    printf("suggestion %s\n", suggestion.c_str());

}

void STD_Suggestion::set_suggestion() {
    if (!SMF["safe"]) {
        suggestion = "restricted";
    } else if (SMF["modern"] && SMF["free"]) {
        suggestion = "strong recommend";
    } else if (SMF["modern"] && !SMF["free"]) {
        suggestion = "recommend";
    } else if (!SMF["modern"] && SMF["free"]) {
        suggestion = "not recommend";
    } else if (!SMF["modern"] && !SMF["free"]) {
        suggestion = "strong not recommend";
    } else {
        std::cerr << "wrong";
    }
    std::cout << std::endl;

}

template<typename T>
struct Thingy {
    T t;
};

Thingy(const char *) -> Thingy<std::string>;


struct factor {
    int f;
};
template<typename term_or_factor_t>
struct term1 {
    term_or_factor_t a;
    factor b;
};

template<typename term_or_factor_t>
term1<term_or_factor_t> operator+(term_or_factor_t a, factor b) {
    return term1<term_or_factor_t>{a, b};
};

template<typename term_t>
int calc(term1<term_t> t) {
    return calc(t.a) + t.b.f;
}

//compling time recursion
template<>
int calc<factor>(term1<factor> t) {
    return t.a.f + t.b.f;
}


//partial specialization
typedef term1<factor> term2_t;
typedef term1<term2_t> term3_t;

template<>
int calc<term2_t>(term3_t t) {
    if (t.b.f == t.a.b.f) {
        std::cout << "transformed" << std::endl;
        return 2 * t.b.f + t.a.a.f;
        //optimization
    } else {
        std::cout << "not translated " << std::endl;
        return t.b.f + t.a.b.f + t.a.a.f;
    }
}

namespace client {
    namespace ast {
        struct nil {
        };
        struct signed_;
        struct program;

        typedef boost::variant<
                nil, unsigned int, boost::recursive_wrapper<signed_>, boost::recursive_wrapper<program>
        >
                operand;

        struct signed_ {
            char sign;
            operand operand_;
        };

        struct operation {
            char operator_;
            operand operand_;
        };

        struct program {
            operand first;
            std::list<operation> rest;
        };
    }
}

BOOST_FUSION_ADAPT_STRUCT(
        client::ast::signed_,
        (char, sign)
                (client::ast::operand, operand_)
)

BOOST_FUSION_ADAPT_STRUCT(
        client::ast::operation,
        (char, operator_)
                (client::ast::operand, operand_)
)

BOOST_FUSION_ADAPT_STRUCT(
        client::ast::program,
        (client::ast::operand, first)
                (std::list<client::ast::operation>, rest)
)

namespace client {
    namespace ast {
        struct printer {
            typedef void result_type;

            void operator()(nil) const {}

            void operator()(unsigned int n) const { std::cout << n; }

            void operator()(operation const &x) const {
                boost::apply_visitor(*this, x.operand_);
                switch (x.operator_) {
                    case '+':
                        std::cout << " add";
                        break;
                    case '-':
                        std::cout << " subt";
                        break;
                    case '*':
                        std::cout << " mult";
                        break;
                    case '/':
                        std::cout << " div";
                        break;
                }
            }

            void operator()(signed_ const &x) const {
                boost::apply_visitor(*this, x.operand_);
                switch (x.sign) {
                    case '-':
                        std::cout << " neg";
                        break;
                    case '+':
                        std::cout << " pos";
                        break;
                }
            }

            void operator()(program const &x) const {
                boost::apply_visitor(*this, x.first);
                BOOST_FOREACH(operation const &oper, x.rest) {
                                std::cout << ' ';
                                (*this)(oper);
                            }
            }
        };

        struct eval {
            typedef int result_type;

            int operator()(nil) const {
                BOOST_ASSERT(0);
                return 0;
            }

            int operator()(unsigned int n) const { return n; }

            int operator()(operation const &x, int lhs) const {
                int rhs = boost::apply_visitor(*this, x.operand_);
                switch (x.operator_) {
                    case '+':
                        return lhs + rhs;
                    case '-':
                        return lhs - rhs;
                    case '*':
                        return lhs * rhs;
                    case '/':
                        return lhs / rhs;
                }
                BOOST_ASSERT(0);
                return 0;
            }

            int operator()(signed_ const &x) const {
                int rhs = boost::apply_visitor(*this, x.operand_);
                switch (x.sign) {
                    case '-':
                        return -rhs;
                    case '+':
                        return +rhs;
                }
                BOOST_ASSERT(0);
                return 0;
            }

            int operator()(program const &x) const {
                int state = boost::apply_visitor(*this, x.first);
                BOOST_FOREACH(operation const &oper, x.rest) {
                                state = (*this)(oper, state);
                            }
                return state;
            }
        };
    }
}

namespace client {
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    template<typename Iterator>
    struct calculator : qi::grammar<Iterator, ast::program(), ascii::space_type> {
        calculator() : calculator::base_type(expression) {
            qi::uint_type uint_;
            qi::char_type char_;

            expression =
                    term
                            >> *((char_('+') >> term)
                                 | (char_('-') >> term)
                            );

            term =
                    factor
                            >> *((char_('*') >> factor)
                                 | (char_('/') >> factor)
                            );

            factor =
                    uint_
                    | '(' >> expression >> ')'
                    | (char_('-') >> factor)
                    | (char_('+') >> factor);
        }

        qi::rule<Iterator, ast::program(), ascii::space_type> expression;
        qi::rule<Iterator, ast::program(), ascii::space_type> term;
        qi::rule<Iterator, ast::operand(), ascii::space_type> factor;
    };
}

int main() {

    STD_Suggestion std1{false, true, true};
    std1.print_suggestion();
    STD_Suggestion std2{true, true, false};
    std2.print_suggestion();

    auto v = std::vector{4, 1, 2, 5};
    for (auto &stda:v) {
        std::cout << stda << std::endl;
    }

    Thingy thing{"A string"};

    factor a{1}, b{2};
    auto t = a + b + b + a;
    // b+b substitute b*2
    auto t2 = a + b + a;

    std::cout << calc(t) << std::endl;

    std::cout << "Pascal's triangle\n";

    for (int n; n < 10; n++) {
        std::cout << std::string(20 - n * 2, ' ');
        for (int k = 1; k < n; ++k) {
            std::cout << std::setw(3) << binom(n, k) << ' ';
        }
        std::cout << '\n';
    }

    typedef client::calculator<std::string::const_iterator> calculator;
    typedef client::ast::program ast_program;
    typedef client::ast::printer ast_print;
    typedef client::ast::eval ast_eval;
    std::string str;

    // calculator interpreter
    while (std::getline(std::cin, str)) {
        if (str.empty() || str[0] == 'q' || str[0] == '0') { break; }

        calculator calc;
        ast_program program;
        ast_print print;
        ast_eval eval;

        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();
        boost::spirit::ascii::space_type space;
        bool r = phrase_parse(iter, end, calc, space, program);

        if (r && iter == end) {
            std::cout << "------------------------------------------\n";
            std::cout << "Parsing succeeded.";
            print(program);
            std::cout << "\nResult: " << eval(program) << std::endl;
            std::cout << "------------------------------------------\n";

        } else {
            std::string rest(iter, end);
            std::cout << "------------------------------------------\n";
            std::cout << "Failed\n";
            std::cout << "Stopped at " << rest << "\n";
            std::cout << "------------------------------------------\n";
        }
    }


    return 0;
}
