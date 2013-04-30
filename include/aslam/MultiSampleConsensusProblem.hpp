#ifndef ASLAM_MULTI_SAMPLE_CONSENSUS_PROBLEM_HPP
#define ASLAM_MULTI_SAMPLE_CONSENSUS_PROBLEM_HPP

#include <boost/random.hpp>
#include <ctime>

namespace aslam {

template<typename MODEL_T>
class MultiSampleConsensusProblem
{
public:
    typedef MODEL_T model_t;
    
    MultiSampleConsensusProblem(bool randomSeed = true);
    virtual ~MultiSampleConsensusProblem();

    virtual void getSamples(int &iterations, std::vector< std::vector<int> > &samples);
    
    virtual bool isSampleGood(const std::vector< std::vector<int> > & sample) const;

    /** \brief Get a pointer to the vector of indices used. */
    boost::shared_ptr< std::vector< std::vector<int> > > getIndices() const;

    void drawIndexSample (std::vector< std::vector<int> > & sample);

    virtual std::vector<int> getSampleSizes() const = 0;

    virtual bool computeModelCoefficients( const std::vector< std::vector<int> > & indices, model_t & outModel) const = 0;

    /** \brief Recompute the model coefficients using the given inlier set
     * and return them to the user. Pure virtual.
     *
     * @note: these are the coefficients of the model after refinement
     * (e.g., after a least-squares optimization)
     *
     * \param[in] inliers the data inliers supporting the model
     * \param[in] model_coefficients the initial guess for the model coefficients
     * \param[out] optimized_coefficients the resultant recomputed coefficients after non-linear optimization
     */
    virtual void optimizeModelCoefficients (const std::vector< std::vector<int> > & inliers,
                                            const model_t & model_coefficients,
                                            model_t & optimized_coefficients) = 0;


    /// \brief evaluate the score for the elements at indices based on this model.
    ///        low scores mean a good fit.
    virtual void getSelectedDistancesToModel( const model_t & model, 
                                              const std::vector< std::vector<int> > & indices,
                                              std::vector< std::vector<double> > & scores) const = 0;



    /** \brief Compute all distances from the cloud data to a given model. Pure virtual.
     * 
     * \param[in] model_coefficients the coefficients of a model that we need to compute distances to 
     * \param[out] distances the resultant estimated distances
     */
    virtual void getDistancesToModel (const model_t & model_coefficients, 
                                      std::vector< std::vector<double> > &distances);



      /** \brief Select all the points which respect the given model
        * coefficients as inliers. Pure virtual.
        * 
        * \param[in] model_coefficients the coefficients of a model that we need to compute distances to
        * \param[in] threshold a maximum admissible distance threshold for determining the inliers from 
        * the outliers
        * \param[out] inliers the resultant model inliers
        */
    virtual void selectWithinDistance (const model_t &model_coefficients, 
                                       const double threshold,
                                       std::vector< std::vector<int> > &inliers);

      /** \brief Count all the points which respect the given model
        * coefficients as inliers. Pure virtual.
        * 
        * \param[in] model_coefficients the coefficients of a model that we need to
        * compute distances to
        * \param[in] threshold a maximum admissible distance threshold for
        * determining the inliers from the outliers
        * \return the resultant number of inliers
        */
      virtual int countWithinDistance (const model_t &model_coefficients, 
                                       const double threshold);



    void setIndices(const std::vector< std::vector<int> > & indices);
    void setUniformIndices(std::vector<int> N);


    int rnd();
    
    int max_sample_checks_;
    
    boost::shared_ptr< std::vector< std::vector<int> > > indices_;
    std::vector< std::vector<int> > shuffled_indices_;

    /** \brief Boost-based random number generator algorithm. */
    boost::mt19937 rng_alg_;

    /** \brief Boost-based random number generator distribution. */
    boost::shared_ptr< boost::uniform_int<> > rng_dist_;

    /** \brief Boost-based random number generator. */
    boost::shared_ptr<boost::variate_generator< boost::mt19937&, boost::uniform_int<> > > rng_gen_;

};

} // namespace aslam

#include "implementation/MultiSampleConsensusProblem.hpp"

#endif /* ASLAM_SAMPLE_CONSENSUS_PROBLEM_HPP */
