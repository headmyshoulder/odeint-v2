#ifndef POLYNOMIAL_HPP_INCLUDED
#define POLYNOMIAL_HPP_INCLUDED

#include <iostream>
#include <boost/array.hpp>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail{

template<
size_t order,
class Time
>
class Polynomial
{
public:
    typedef Time time_type;
    typedef Polynomial<order, Time> poly_type;

    Polynomial()
    {
        for(size_t i = 0; i<order; ++i)
            m_coeff[i] = 0;
        // start with 'constant' polynomial
        m_coeff[order-1] = 1;
    };

    Polynomial(boost::array<time_type, order> coeff)
    {
        m_coeff = coeff;
    };

    Polynomial<order+1, time_type> integrate()
    {
        boost::array<time_type, order + 1> coeff_int;
        coeff_int[order] = 0;

        for(size_t i=0; i<order; ++i)
        {
            if(i != order)
                coeff_int[i] = m_coeff[i] / (order - i);
        }

        Polynomial<order+1, time_type> polyInt(coeff_int);
        return polyInt;
    };

    time_type evaluate(const time_type &t)
    {
        // fma: x*y+z
        m_res = m_coeff[0];
        for(size_t i=1; i<order; ++i)
        {
            // m_res = fma(m_res, t, m_coeff[i]);
            m_res = m_coeff[i] + m_res * t;
        }

        return m_res;
    };

    time_type evaluate_integrated(const time_type &t)
    {
        m_res = m_coeff[0]/order;
        for(size_t i=1; i<order+1; ++i)
        {
            // m_res = fma(m_res, t, ((i >= order)?0:m_coeff[i]/(order-i)));
            m_res = ((i >= order)?0:m_coeff[i]/(order-i)) + m_res * t;
        }

        return m_res;
    }

    void add_root(const time_type &root)
    {
        for(size_t j=0; j<order; ++j) // updating all coefficients
        {
            m_coeff[j] = -m_coeff[j]*root + ((j<order-1)?m_coeff[j+1]:0);
        }
    };

    void reset()
    {
        for(size_t i = 0; i<order; ++i)
            m_coeff[i] = 0;
        
        m_coeff[order-1] = 1;
    };

    void remove_root(const time_type &root)
    {
        // 'reverse' add_root; synthetic division
        time_type tmp[2];

        m_coeff[0] = 0;
        for(size_t j=0; j<order; ++j)
        {
            tmp[0] = (j!=0)?tmp[1]:0;
            tmp[1] = m_coeff[j+1];

            m_coeff[j+1] = tmp[0] + root * m_coeff[j];
        }
    };

    void pretty_print()
    {
        for(size_t i=0; i<order; ++i)
        {
            std::cout << m_coeff[i];
            std::cout << "x^" << order - i - 1;
            if(i != order - 1)
                std::cout<< " + ";
            else
                std::cout << std::endl;
        }
    };

    static Polynomial<order, time_type> from_roots(time_type * roots, unsigned short num_roots)
    {
        time_type coeff[order] = {0};
        coeff[order-1] = 1;

        for(size_t i=0; i<num_roots; ++i) // going over all roots
        {
            for(size_t j=0; j<order; ++j) // updating all coefficients
            {
                coeff[j] = -coeff[j]*roots[i] + ((j<order)?coeff[j+1]:0);
            }
        }
        Polynomial<order, time_type> poly(coeff);
        return poly;
    };

    // first element is highest order
    boost::array<time_type, order> m_coeff;
private:
    time_type m_res;
};

}
}
}
}

#endif