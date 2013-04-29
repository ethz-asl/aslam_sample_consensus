#ifndef ASLAM_MULTI_RANSAC_HPP
#define ASLAM_MULTI_RANSAC_HPP

#include <boost/shared_ptr.hpp>
#include <vector>
#include <aslam/MultiSampleConsensus.hpp>
#include <cstdio>

namespace aslam {
    
    /// \brief plain old ransac
    template<typename PROBLEM_T>
    class MultiRansac : public MultiSampleConsensus<PROBLEM_T>
    {
    public:
        typedef PROBLEM_T problem_t;
        typedef typename problem_t::model_t model_t;
        using MultiSampleConsensus<problem_t>::max_iterations_;
        using MultiSampleConsensus<problem_t>::threshold_;
        using MultiSampleConsensus<problem_t>::iterations_;
        using MultiSampleConsensus<problem_t>::sac_model_;
        using MultiSampleConsensus<problem_t>::model_;
        using MultiSampleConsensus<problem_t>::model_coefficients_;
        using MultiSampleConsensus<problem_t>::inliers_;
        using MultiSampleConsensus<problem_t>::probability_;
        
        MultiRansac(int maxIterations = 1000, double threshold = 1.0, double probability = 0.99);
        virtual ~MultiRansac();
        
        bool computeModel(int debug_verbosity_level = 0);

    };

} // namespace aslam

#include "implementation/MultiRansac.hpp"

#endif /* ASLAM_MULTI_RANSAC_HPP */
