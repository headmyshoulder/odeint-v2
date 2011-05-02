#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/sequence/intrinsic/end.hpp>
#include <boost/fusion/iterator/equal_to.hpp>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/iterator/deref.hpp>
#include <boost/fusion/iterator/distance.hpp>
#include <boost/mpl/bool.hpp>

namespace boost { namespace fusion {
namespace detail
{
	template<>
    struct for_each_unrolled<6>
    {
        template<typename I0, typename F>
        static void call(I0 const& i0, F const& f)
        {
            f(*i0);
            typedef typename result_of::next<I0>::type I1;
            I1 i1(fusion::next(i0));
            f(*i1);
            typedef typename result_of::next<I1>::type I2;
            I2 i2(fusion::next(i1));
            f(*i2);
			typedef typename result_of::next<I2>::type I3;
            I3 i3(fusion::next(i2));
            f(*i3);
			typedef typename result_of::next<I3>::type I4;
            I4 i4(fusion::next(i3));
            f(*i4);
			typedef typename result_of::next<I4>::type I5;
            I5 i5(fusion::next(i4));
            f(*i5);
        }
    };
}
} }