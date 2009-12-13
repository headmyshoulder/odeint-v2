/* Boost odeint/detail/accumulators.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner
 
 Some algebraic operations for iterators

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP
#define BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP

#include <iostream>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {
namespace it_algebra { // iterator algebra


    // computes y += alpha * x1
    template <
        class InOutIterator ,
        class InputIterator ,
        class T
        >
    void increment(
                   InOutIterator first1 ,
                   InOutIterator last1 ,
                   InputIterator first2 ,
                   T alpha
                   )
    {
        while( first1 != last1 )
            (*first1++) += alpha * (*first2++);
    }

    // computes y = x1 - x2
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2
        >
    void assign_diff(
                     OutputIterator first1 ,
                     OutputIterator last1 ,
                     InputIterator1 first2 ,
                     InputIterator2 first3 )
    {
        while( first1 != last1 )
            (*first1++) = (*first2++) - (*first3++);
    }


    // computes y = x1 + alpha * x2
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class T
        >
    void assign_sum(
                    OutputIterator first1 ,
                    OutputIterator last1 ,
                    InputIterator1 first2 ,
                    InputIterator2 first3 ,
                    T alpha )
    {
        while( first1 != last1 )
            (*first1++) = (*first2++) + alpha * (*first3++);
    }



    // computes y = alpha1 * ( x1 + x2 + alpha2*x3 )
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class T
        >
    void increment_sum_sum(
                           OutputIterator first1 ,
                           OutputIterator last1 ,
                           InputIterator1 first2 ,
                           InputIterator2 first3 ,
                           InputIterator3 first4 ,
                           T alpha1 ,
                           T alpha2
                           )
    {
        while( first1 != last1 )
            (*first1++) += alpha1 *
                ( (*first2++) + (*first3++) + alpha2 * (*first4++) );
    }




    // computes y = alpha1*x1 + alpha2*x2
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin )
    {
        while( y_begin != y_end )
        {
            (*y_begin++) = alpha1 * (*x1_begin++) +
                alpha2 * (*x2_begin++);
        }
    }


    // computes y = x1 + alpha2*x2 + alpha3*x3
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin )
    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++);
    }

    // computes y = x1 + alpha2*x2 + alpha3*x3 + alpha4*x4
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class InputIterator4 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin ,
                           T alpha4 ,
                           InputIterator4 x4_begin )
    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++);
    }

    // computes y = x1 + alpha2*x2 + alpha3*x3 + alpha4*x4 + alpha5*x5
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class InputIterator4 ,
        class InputIterator5 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin ,
                           T alpha4 ,
                           InputIterator4 x4_begin ,
                           T alpha5 ,
                           InputIterator5 x5_begin )
    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++) +
                alpha5 * (*x5_begin++);
    }


    // computes y = x1 + alpha2*x2 + alpha3*x3 + alpha4*x4 + alpha5*x5
    // + alpha6*x6
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class InputIterator4 ,
        class InputIterator5 ,
        class InputIterator6 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin ,
                           T alpha4 ,
                           InputIterator4 x4_begin ,
                           T alpha5 ,
                           InputIterator5 x5_begin ,
                           T alpha6 ,
                           InputIterator6 x6_begin )
    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++) +
                alpha5 * (*x5_begin++) +
                alpha6 * (*x6_begin++);
    }


    // computes y = x1 + alpha2*x2 + alpha3*x3 + alpha4*x4 + alpha5*x5
    // + alpha6*x6 + alpha7*x7
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class InputIterator4 ,
        class InputIterator5 ,
        class InputIterator6 ,
        class InputIterator7 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin ,
                           T alpha4 ,
                           InputIterator4 x4_begin ,
                           T alpha5 ,
                           InputIterator5 x5_begin ,
                           T alpha6 ,
                           InputIterator6 x6_begin ,
                           T alpha7 ,
                           InputIterator7 x7_begin )

    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++) +
                alpha5 * (*x5_begin++) +
                alpha6 * (*x6_begin++) +
                alpha7 * (*x7_begin++) ;
    }



    // computes y = x1 + alpha2*x2 + alpha3*x3 + alpha4*x4 + alpha5*x5
    // + alpha6*x6 + alpha7*x7 + alpha8*x8
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class InputIterator4 ,
        class InputIterator5 ,
        class InputIterator6 ,
        class InputIterator7 ,
        class InputIterator8 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin ,
                           T alpha4 ,
                           InputIterator4 x4_begin ,
                           T alpha5 ,
                           InputIterator5 x5_begin ,
                           T alpha6 ,
                           InputIterator6 x6_begin ,
                           T alpha7 ,
                           InputIterator7 x7_begin ,
                           T alpha8 ,
                           InputIterator8 x8_begin )
    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++) +
                alpha5 * (*x5_begin++) +
                alpha6 * (*x6_begin++) +
                alpha7 * (*x7_begin++) +
                alpha8 * (*x8_begin++);
    }

    // computes y = x1 + alpha2*x2 + alpha3*x3 + alpha4*x4 + alpha5*x5
    // + alpha6*x6 + alpha7*x7 + alpha8*x8 + alpha9*x9
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,        
        class InputIterator4 ,
        class InputIterator5 ,
        class InputIterator6 ,
        class InputIterator7 ,
        class InputIterator8 ,
        class InputIterator9 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin ,
                           T alpha4 ,
                           InputIterator4 x4_begin ,
                           T alpha5 ,
                           InputIterator5 x5_begin ,
                           T alpha6 ,
                           InputIterator6 x6_begin ,
                           T alpha7 ,
                           InputIterator7 x7_begin ,
                           T alpha8 ,
                           InputIterator8 x8_begin ,
                           T alpha9 ,
                           InputIterator9 x9_begin )
    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++) +
                alpha5 * (*x5_begin++) +
                alpha6 * (*x6_begin++) +
                alpha7 * (*x7_begin++) +
                alpha8 * (*x8_begin++) +
                alpha9 * (*x9_begin++);
    }




    // computes y = x1 + alpha2*x2 + alpha3*x3 + alpha4*x4 + alpha5*x5
    // + alpha6*x6 + alpha7*x7 + alpha8*x8 + alpha9*x9 + alpha10*x10
    template <
        class OutputIterator ,
        class InputIterator1 ,
        class InputIterator2 ,
        class InputIterator3 ,
        class InputIterator4 ,
        class InputIterator5 ,
        class InputIterator6 ,
        class InputIterator7 ,
        class InputIterator8 ,
        class InputIterator9 ,
        class InputIterator10 ,
        class T
        >
    inline void scale_sum( OutputIterator y_begin ,
                           OutputIterator y_end ,
                           T alpha1 ,
                           InputIterator1 x1_begin ,
                           T alpha2 ,
                           InputIterator2 x2_begin ,
                           T alpha3 ,
                           InputIterator3 x3_begin ,
                           T alpha4 ,
                           InputIterator4 x4_begin ,
                           T alpha5 ,
                           InputIterator5 x5_begin ,
                           T alpha6 ,
                           InputIterator6 x6_begin ,
                           T alpha7 ,
                           InputIterator7 x7_begin ,
                           T alpha8 ,
                           InputIterator8 x8_begin ,
                           T alpha9 ,
                           InputIterator9 x9_begin ,
                           T alpha10 ,
                           InputIterator10 x10_begin 
        )
    {
        while( y_begin != y_end )
            (*y_begin++) = 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++) +
                alpha5 * (*x5_begin++) +
                alpha6 * (*x6_begin++) +
                alpha7 * (*x7_begin++) +
                alpha8 * (*x8_begin++) +
                alpha9 * (*x9_begin++) +
                alpha10 * (*x10_begin++);
    }







    // generic version for n values
    template <
        class OutputIterator ,
        class InputIterator ,
        class InputIteratorIterator ,
        class FactorIterator ,
        class T
        >
    inline void scale_sum_generic( OutputIterator y_begin ,
                                   OutputIterator y_end ,
                                   FactorIterator alpha_begin ,
                                   FactorIterator alpha_end ,
                                   T beta ,
                                   InputIterator x_begin ,
                                   InputIteratorIterator x_iter_begin )
    {
        FactorIterator alpha_iter;
        InputIteratorIterator x_iter_iter;
        while( y_begin != y_end )
	{
            x_iter_iter = x_iter_begin;
            alpha_iter = alpha_begin;
            *y_begin = *x_begin++;
            //std::clog<<(*y_begin);
            while( alpha_iter != alpha_end )
	    {
                //std::clog<< " + " <<beta<<" * "<<*alpha_iter<<"*"<<(*(*(x_iter_iter)));
                (*y_begin) += beta * (*alpha_iter++) * (*(*x_iter_iter++)++);
            }
            //std::clog<<" = "<<*y_begin<<std::endl;
            y_begin++;
        }
        //std::clog<<std::endl;
    }


    // computes y += alpha1*x1 + alpha2*x2 + alpha3*x3 + alpha4*x4
    template <
        class OutputIterator,
        class InputIterator1,
        class InputIterator2,
        class InputIterator3,
        class InputIterator4,
        class T
        >
    inline void scale_sum_inplace(
            OutputIterator y_begin,
            OutputIterator y_end,
            T alpha1,
            InputIterator1 x1_begin,
            T alpha2,
            InputIterator2 x2_begin,
            T alpha3,
            InputIterator3 x3_begin,
            T alpha4,
            InputIterator4 x4_begin )
    {   
        while( y_begin != y_end )
            (*y_begin++) += 
                alpha1 * (*x1_begin++) + 
                alpha2 * (*x2_begin++) + 
                alpha3 * (*x3_begin++) +
                alpha4 * (*x4_begin++);
    }

    /* calculates tmp = y, y = x1 + alpha*x2, x1 = tmp */
    template <
        class OutputIterator,
        class InputIterator,
        class T
        >
    inline void scale_sum_swap(
            OutputIterator y_begin,
            OutputIterator y_end,
            OutputIterator x1_begin,
            T alpha,
            InputIterator x2_begin )
    {
        T swap;
        while( y_begin != y_end )
        {
            swap = (*x1_begin) + alpha*(*x2_begin++);
            *x1_begin++ = *y_begin;
            *y_begin++ = swap;
        }
    }


    // computes y = x1 + alpha2 * x2 ; x2 += x3
    template<
        class OutputIterator ,
        class InputIterator1 ,
        class InOutIterator ,
        class InputIterator2 ,
        class T
        >
    void assign_sum_increment(
                              OutputIterator first1 ,
                              OutputIterator last1 ,
                              InputIterator1 first2 ,
                              InOutIterator first3 ,
                              InputIterator2 first4 ,
                              T alpha
                              )
    {
        while( first1 != last1 )
          {
                (*first1++) = (*first2++) + alpha * (*first3);
                (*first3++) += (*first4++);
          }
    }

    
    




}
}
}
}
}


#endif //BOOST_NUMERIC_ODEINT_DETAIL_ACCUMULATORS_HPP
