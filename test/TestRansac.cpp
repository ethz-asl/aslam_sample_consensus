#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <gtest/gtest.h>

#include <aslam/Ransac.hpp>
#include <aslam/Msac.hpp>
//#include <aslam/Rmsac.hpp>
#include <aslam/Lmeds.hpp>
#include <aslam/Prosac.hpp>
#include <Eigen/Dense>
#include <aslam/SampleConsensusProblem.hpp>





class CircleFitProblem : public aslam::SampleConsensusProblem<Eigen::VectorXd>
{
public:
    typedef Eigen::VectorXd model_t;

    CircleFitProblem(){}
    virtual ~CircleFitProblem(){}

    /// \brief how many elements in the problem.
    size_t numElements() const{ return points_.cols(); }
        
    virtual int getSampleSize() const { return 3; }
        

    virtual bool computeModelCoefficients( const std::vector<int> & samples, model_t & model_coefficients) const
        {
            // Need 3 samples
            if (samples.size () != 3)
            {
                fprintf(stderr,"[sm::SampleConsensusModelCircle2D::computeModelCoefficients] Invalid set of samples given (%zu)!\n", samples.size ());
                return (false);
            }

            model_coefficients.resize (3);

            Eigen::Vector2d p0 (points_(0,samples[0]), points_(1,samples[0]));
            Eigen::Vector2d p1 (points_(0,samples[1]), points_(1,samples[1]));
            Eigen::Vector2d p2 (points_(0,samples[2]), points_(1,samples[2]));

            Eigen::Vector2d u = (p0 + p1) / 2.0;
            Eigen::Vector2d v = (p1 + p2) / 2.0;

            Eigen::Vector2d p1p0dif = p1 - p0;
            Eigen::Vector2d p2p1dif = p2 - p1;
            Eigen::Vector2d uvdif   = u - v;

            Eigen::Vector2d m (- p1p0dif[0] / p1p0dif[1], - p2p1dif[0] / p2p1dif[1]);

            // Center (x, y)
            model_coefficients[0] = static_cast<float> ((m[0] * u[0] -  m[1] * v[0]  - uvdif[1] )             / (m[0] - m[1]));
            model_coefficients[1] = static_cast<float> ((m[0] * m[1] * uvdif[0] +  m[0] * v[1] - m[1] * u[1]) / (m[0] - m[1]));

            // Radius
            model_coefficients[2] = static_cast<float> (sqrt ((model_coefficients[0] - p0[0]) * (model_coefficients[0] - p0[0]) +
                                                              (model_coefficients[1] - p0[1]) * (model_coefficients[1] - p0[1])));
            return (true);
        }


    virtual void optimizeModelCoefficients (const std::vector<int> & inliers, 
                                            const model_t & model_coefficients,
                                            model_t & optimized_coefficients)
        {
            optimized_coefficients = model_coefficients;
        }


    /// \brief evaluate the score for the elements at indices based on this model.
    ///        low scores mean a good fit.
    virtual void getSelectedDistancesToModel( const model_t & model_coefficients, 
                                              const std::vector<int> & indices,
                                              std::vector<double> & scores) const
        {
            scores.resize(indices.size());
            // Iterate through the 3d points and calculate the distances from them to the sphere
            for (size_t i = 0; i < indices.size (); ++i)
            {
                // Calculate the distance from the point to the sphere as the difference between
                // dist(point,sphere_origin) and sphere_radius
                scores[i] = fabsf (sqrtf (
                                       ( points_(0,indices[i]) - model_coefficients[0] ) *
                                       ( points_(0,indices[i]) - model_coefficients[0] ) +
                                       
                                       ( points_(1,indices[i]) - model_coefficients[1] ) *
                                       ( points_(1,indices[i]) - model_coefficients[1] )
                                       ) - model_coefficients[2]);
            }
        }
    
    
    Eigen::MatrixXd points_;

        
};


boost::shared_ptr<CircleFitProblem> genCfp()
{
    boost::shared_ptr<CircleFitProblem> cfp_ptr( new CircleFitProblem );
    CircleFitProblem & cfp = *cfp_ptr;

    const int N = 100;
    const int G = 80;

    Eigen::Vector2d c(5,10);
    double r = 2;

    //sm::console::setVerbosityLevel (sm::console::L_DEBUG);


    cfp.points_.resize(2,N);
    for(int i = 0; i < N; ++i)
    {
        Eigen::Vector2d p;
        p.setRandom();
        p[0] -= 0.5;
        p[1] -= 0.5;
  
        if(i < G)
        {
            p.normalize();
            cfp.points_.col(i) = p * r + c;
        }
        else
        {
            cfp.points_.col(i) = p * r;
        }
    }
    
    cfp.setUniformIndices(N);

    return cfp_ptr;
}


TEST(RansacTestSuite, testRansac)
{
    boost::shared_ptr<CircleFitProblem> cfp_ptr = genCfp();
    aslam::Ransac<CircleFitProblem> ransac(100, 0.05, 0.99);
    
    ransac.sac_model_ = cfp_ptr;
    ransac.computeModel(4);

    fprintf(stdout, "Model %f, %f, %f\n", ransac.model_coefficients_[0], ransac.model_coefficients_[1], ransac.model_coefficients_[2]);
    fprintf(stdout, "Iterations: %d\n", ransac.iterations_);    
    fprintf(stdout, "Inliers: %zu\n", ransac.inliers_.size());

    ASSERT_NEAR(ransac.model_coefficients_[0], 5.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[1], 10.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[2], 2.0, 1e-2);

}


TEST(RansacTestSuite, testMsac)
{
    boost::shared_ptr<CircleFitProblem> cfp_ptr = genCfp();
    aslam::Msac<CircleFitProblem> ransac(100, 0.05, 0.99);
    
    ransac.sac_model_ = cfp_ptr;
    ransac.computeModel(4);

    fprintf(stdout, "Model %f, %f, %f\n", ransac.model_coefficients_[0], ransac.model_coefficients_[1], ransac.model_coefficients_[2]);
    fprintf(stdout, "Iterations: %d\n", ransac.iterations_);    
    fprintf(stdout, "Inliers: %zu\n", ransac.inliers_.size());

    ASSERT_NEAR(ransac.model_coefficients_[0], 5.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[1], 10.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[2], 2.0, 1e-2);

}

TEST(RansacTestSuite, testLmeds)
{
    boost::shared_ptr<CircleFitProblem> cfp_ptr = genCfp();
    aslam::Lmeds<CircleFitProblem> ransac(100, 0.05, 0.99);
    
    ransac.sac_model_ = cfp_ptr;
    ransac.computeModel(4);

    fprintf(stdout, "Model %f, %f, %f\n", ransac.model_coefficients_[0], ransac.model_coefficients_[1], ransac.model_coefficients_[2]);
    fprintf(stdout, "Iterations: %d\n", ransac.iterations_);    
    fprintf(stdout, "Inliers: %zu\n", ransac.inliers_.size());

    ASSERT_NEAR(ransac.model_coefficients_[0], 5.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[1], 10.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[2], 2.0, 1e-2);

}


TEST(RansacTestSuite, testProsac)
{
    boost::shared_ptr<CircleFitProblem> cfp_ptr = genCfp();
    aslam::Prosac<CircleFitProblem> ransac(100, 0.05, 0.99);
    
    ransac.sac_model_ = cfp_ptr;
    ransac.computeModel(4);

    fprintf(stdout, "Model %f, %f, %f\n", ransac.model_coefficients_[0], ransac.model_coefficients_[1], ransac.model_coefficients_[2]);
    fprintf(stdout, "Iterations: %d\n", ransac.iterations_);    
    fprintf(stdout, "Inliers: %zu\n", ransac.inliers_.size());

    ASSERT_NEAR(ransac.model_coefficients_[0], 5.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[1], 10.0, 1e-2);
    ASSERT_NEAR(ransac.model_coefficients_[2], 2.0, 1e-2);

}


// TEST(RansacTestSuite, testRmsac)
// {
//     boost::shared_ptr<CircleFitProblem> cfp_ptr = genCfp();
//     aslam::Rmsac<CircleFitProblem> ransac(100, 0.05, 0.99);
    
//     ransac.sac_model_ = cfp_ptr;
//     ransac.computeModel(4);

//     fprintf(stdout, "Model %f, %f, %f\n", ransac.model_coefficients_[0], ransac.model_coefficients_[1], ransac.model_coefficients_[2]);
//     fprintf(stdout, "Iterations: %d\n", ransac.iterations_);    
//     fprintf(stdout, "Inliers: %zu\n", ransac.inliers_.size());

//     ASSERT_NEAR(ransac.model_coefficients_[0], 5.0, 1e-2);
//     ASSERT_NEAR(ransac.model_coefficients_[1], 10.0, 1e-2);
//     ASSERT_NEAR(ransac.model_coefficients_[2], 2.0, 1e-2);

// }

