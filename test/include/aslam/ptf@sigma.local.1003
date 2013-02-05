#ifndef ASLAM_MLESAC_HPP
#define ASLAM_MLESAC_HPP

#include <sm/console.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <aslam/SampleConsensus.hpp>

namespace aslam {
    
    /// \brief plain old mlesac
    template<typename PROBLEM_T>
    class Mlesac : public SampleConsensus<PROBLEM_T>
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
        
        Mlesac(int maxIterations = 1000, double threshold = 1.0, double probability = 0.99);
        virtual ~Mlesac();
        
        bool computeModel(int debug_verbosity_level = 0);
        /** \brief Maximum number of EM (Expectation Maximization) iterations. */
        int iterations_EM_;
        /** \brief The MLESAC sigma parameter. */
        double sigma_;
    };

} // namespace aslam

#include "implementation/Mlesac.hpp"

#endif /* ASLAM_MLESAC_HPP */
