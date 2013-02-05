#ifndef ASLAM_PROSAC_HPP
#define ASLAM_PROSAC_HPP

#include <boost/shared_ptr.hpp>
#include <vector>
#include <aslam/SampleConsensus.hpp>
#include <cmath>
#include <limits>
#include <boost/math/distributions/binomial.hpp>

namespace aslam {
    
    /// \brief plain old prosac
    template<typename PROBLEM_T>
    class Prosac : public SampleConsensus<PROBLEM_T>
    {
    public:
        typedef PROBLEM_T problem_t;
        typedef typename problem_t::model_t model_t;
        using SampleConsensus<problem_t>::max_iterations_;
        using SampleConsensus<problem_t>::threshold_;
        using SampleConsensus<problem_t>::iterations_;
        using SampleConsensus<problem_t>::sac_model_;
        using SampleConsensus<problem_t>::model_;
        using SampleConsensus<problem_t>::model_coefficients_;
        using SampleConsensus<problem_t>::inliers_;
        using SampleConsensus<problem_t>::probability_;
        
        Prosac(int maxIterations = 1000, double threshold = 1.0, double probability = 0.99);
        virtual ~Prosac();
        
        bool computeModel(int debug_verbosity_level = 0);

    };

} // namespace aslam

#include "implementation/Prosac.hpp"

#endif /* ASLAM_PROSAC_HPP */
