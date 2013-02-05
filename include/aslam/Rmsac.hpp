#ifndef ASLAM_RMSAC_HPP
#define ASLAM_RMSAC_HPP

#include <sm/console.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <aslam/SampleConsensus.hpp>
#include <sm/round.hpp>

namespace aslam {
    
    /// \brief plain old rmsac
    template<typename PROBLEM_T>
    class Rmsac : public SampleConsensus<PROBLEM_T>
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
        
        Rmsac(int maxIterations = 1000, double threshold = 1.0, double probability = 0.99);
        virtual ~Rmsac();
        
        bool computeModel(int debug_verbosity_level = 0);

      /** \brief Set the percentage of points to pre-test.
        * \param nr_pretest percentage of points to pre-test
        */
      inline void setFractionNrPretest (double nr_pretest) { fraction_nr_pretest_ = nr_pretest; }

      /** \brief Get the percentage of points to pre-test. */
      inline double getFractionNrPretest () { return (fraction_nr_pretest_); }


        /** \brief Number of samples to randomly pre-test, in percents. */
        double fraction_nr_pretest_;

    };

} // namespace aslam

#include "implementation/Rmsac.hpp"

#endif /* ASLAM_RMSAC_HPP */
