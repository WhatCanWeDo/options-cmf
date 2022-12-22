#include <boost/math/distributions/normal.hpp>
#include <cmath>

using boost::math::normal;

struct stock {
    double price;
    double vol;
    double dividend_rate; // annualized

    stock() = delete;

    stock(double price, double vol) 
        : price(price), vol(vol), dividend_rate(0) {}

    stock(double price, double vol, double dividend_rate)
        : price(price), vol(vol), dividend_rate(dividend_rate) {}
};

struct option {
    double strike;
    double ttm;
    stock st;
};

struct model {
    double risk_free_rate; // annualized
    option opt;

    model(double risk_free_rate, option opt)
        : risk_free_rate(risk_free_rate), opt(std::move(opt)) {
            d1 = (log(opt.st.price / opt.strike) \
                + (risk_free_rate - opt.st.dividend_rate + (opt.st.vol * opt.st.vol) / 2) * opt.ttm) \
                / (opt.st.vol * sqrt(opt.ttm));
            d2 = d1 - sqrt(opt.ttm) * opt.st.vol;
        }

    double price_call() {
        return opt.st.price * cdf(dist, d1) * exp(-opt.st.dividend_rate * opt.ttm) \
                - opt.strike * cdf(dist, d2) * exp(-risk_free_rate * opt.ttm);
    }

    double price_put() {
        return opt.strike * cdf(dist, -d2) * exp(-risk_free_rate * opt.ttm) \
                - opt.st.price * cdf(dist, -d1) * exp(-opt.st.dividend_rate * opt.ttm);
    }

    double delta() {
        return exp(-opt.st.dividend_rate * opt.ttm) * cdf(dist, -d1);
    }

    double vega() {
        return opt.st.price * exp(-opt.st.dividend_rate * opt.ttm) * sqrt(opt.ttm) * pdf(dist, d1);
    }

    double psi() {
        return -opt.st.price * exp(-opt.st.dividend_rate * opt.ttm) \
            * (sqrt(opt.ttm)/opt.st.vol * pdf(dist, d1) + opt.ttm * cdf(dist, -d1)) \
            + opt.strike * exp(-risk_free_rate * opt.ttm) * pdf(dist, -d2) \
            * sqrt(opt.ttm) / opt.st.vol;
    }

    double theta() {
        return -opt.st.price * pdf(dist, d1) * opt.st.vol / (2 * sqrt(opt.ttm)) \
            - opt.strike * risk_free_rate * exp(-risk_free_rate * opt.ttm) * cdf(dist, -d2) \
            + opt.st.dividend_rate * opt.st.price * exp(-opt.st.dividend_rate * opt.ttm) * cdf(dist, -d1);
    }

    double rho() {
        return -opt.strike * opt.ttm * exp(-risk_free_rate * opt.ttm) * cdf(dist, -d2);
    }

    double gamma() {
        return exp(-opt.st.dividend_rate * opt.ttm) * pdf(dist, d1) / (opt.st.price * opt.st.vol * sqrt(opt.ttm));
    }

    double volga() {
        return vega() * d1 * d2 / opt.st.vol;
    }
private:
    double d1, d2;
    normal dist;
};
