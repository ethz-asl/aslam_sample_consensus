
namespace aslam {

    template<typename P>
    Prosac<P>::Prosac(int maxIterations, double threshold, double probability) :
        SampleConsensus<P>(maxIterations, threshold, probability)
    {

    }

    template<typename P>
    Prosac<P>::~Prosac(){}


    template<typename PROBLEM_T>
    bool Prosac<PROBLEM_T>::computeModel(int debug_verbosity_level)
    {
        typedef PROBLEM_T problem_t;
        typedef typename problem_t::model_t model_t;

        // Warn and exit if no threshold was set
        if (threshold_ == std::numeric_limits<double>::max())
        {
            SM_ERROR ("[sm::ProgressiveSampleConsensus::computeModel] No threshold set!\n");
            return (false);
        }

        // Initialize some PROSAC constants
        const int T_N = 200000;
        const size_t N = sac_model_->indices_->size ();
        const size_t m = sac_model_->getSampleSize ();
        float T_n = static_cast<float> (T_N);
        for (unsigned int i = 0; i < m; ++i)
            T_n *= static_cast<float> (m - i) / static_cast<float> (N - i);
        float T_prime_n = 1.0f;
        size_t I_N_best = 0;
        float n = static_cast<float> (m);

        // Define the n_Start coefficients from Section 2.2
        float n_star = static_cast<float> (N);
        float epsilon_n_star = 0.0;
        size_t k_n_star = T_N;

        // Compute the I_n_star_min of Equation 8
        std::vector<unsigned int> I_n_star_min (N);

        // Initialize the usual RANSAC parameters
        iterations_ = 0;

        std::vector<int> inliers;
        std::vector<int> selection;
        model_t model_coefficients;

        // We will increase the pool so the indices_ vector can only contain m elements at first
        std::vector<int> index_pool;
        index_pool.reserve (N);
        for (unsigned int i = 0; i < n; ++i)
            index_pool.push_back (sac_model_->indices_->operator[](i));

        // Iterate
        while (static_cast<unsigned int> (iterations_) < k_n_star)
        {
            // Choose the samples

            // Step 1
            // According to Equation 5 in the text text, not the algorithm
            if ((iterations_ == T_prime_n) && (n < n_star))
            {
                // Increase the pool
                ++n;
                if (n >= N)
                    break;
                index_pool.push_back (sac_model_->indices_->at(static_cast<unsigned int> (n - 1)));
                // Update other variables
                float T_n_minus_1 = T_n;
                T_n *= (static_cast<float>(n) + 1.0f) / (static_cast<float>(n) + 1.0f - static_cast<float>(m));
                T_prime_n += ceilf (T_n - T_n_minus_1);
            }

            // Step 2
            sac_model_->indices_->swap (index_pool);
            selection.clear ();
            sac_model_->getSamples (iterations_, selection);
            if (T_prime_n < iterations_)
            {
                selection.pop_back ();
                selection.push_back (sac_model_->indices_->at(static_cast<unsigned int> (n - 1)));
            }

            // Make sure we use the right indices for testing
            sac_model_->indices_->swap (index_pool);

            if (selection.empty ())
            {
                SM_ERROR ("[sm::ProgressiveSampleConsensus::computeModel] No samples could be selected!\n");
                break;
            }

            // Search for inliers in the point cloud for the current model
            if (!sac_model_->computeModelCoefficients (selection, model_coefficients))
            {
                ++iterations_;
                continue;
            }

            // Select the inliers that are within threshold_ from the model
            inliers.clear ();
            sac_model_->selectWithinDistance (model_coefficients, threshold_, inliers);

            size_t I_N = inliers.size ();

            // If we find more inliers than before
            if (I_N > I_N_best)
            {
                I_N_best = I_N;

                // Save the current model/inlier/coefficients selection as being the best so far
                inliers_ = inliers;
                model_ = selection;
                model_coefficients_ = model_coefficients;

                // We estimate I_n_star for different possible values of n_star by using the inliers
                std::sort (inliers.begin (), inliers.end ());

                // Try to find a better n_star
                // We minimize k_n_star and therefore maximize epsilon_n_star = I_n_star / n_star
                size_t possible_n_star_best = N, I_possible_n_star_best = I_N;
                float epsilon_possible_n_star_best = static_cast<float>(I_possible_n_star_best) / static_cast<float>(possible_n_star_best);

                // We only need to compute possible better epsilon_n_star for when _n is just about to be removed an inlier
                size_t I_possible_n_star = I_N;
                for (std::vector<int>::const_reverse_iterator last_inlier = inliers.rbegin (), 
                         inliers_end = inliers.rend (); 
                     last_inlier != inliers_end; 
                     ++last_inlier, --I_possible_n_star)
                {
                    // The best possible_n_star for a given I_possible_n_star is the index of the last inlier
                    unsigned int possible_n_star = (*last_inlier) + 1;
                    if (possible_n_star <= m)
                        break;

                    // If we find a better epsilon_n_star
                    float epsilon_possible_n_star = static_cast<float>(I_possible_n_star) / static_cast<float>(possible_n_star);
                    // Make sure we have a better epsilon_possible_n_star
                    if ((epsilon_possible_n_star > epsilon_n_star) && (epsilon_possible_n_star > epsilon_possible_n_star_best))
                    {
                        using namespace boost::math;
                        // Typo in Equation 7, not (n-m choose i-m) but (n choose i-m)
                        size_t I_possible_n_star_min = m
                            + static_cast<size_t> (ceil (quantile (complement (binomial_distribution<float>(static_cast<float> (possible_n_star), 0.1f), 0.05))));
                        // If Equation 9 is not verified, exit
                        if (I_possible_n_star < I_possible_n_star_min)
                            break;

                        possible_n_star_best = possible_n_star;
                        I_possible_n_star_best = I_possible_n_star;
                        epsilon_possible_n_star_best = epsilon_possible_n_star;
                    }
                }

                // Check if we get a better epsilon
                if (epsilon_possible_n_star_best > epsilon_n_star)
                {
                    // update the best value
                    epsilon_n_star = epsilon_possible_n_star_best;

                    // Compute the new k_n_star
                    float bottom_log = 1 - std::pow (epsilon_n_star, static_cast<float>(m));
                    if (bottom_log == 0)
                        k_n_star = 1;
                    else if (bottom_log == 1)
                        k_n_star = T_N;
                    else
                        k_n_star = static_cast<int> (ceil (log (0.05) / log (bottom_log)));
                    // It seems weird to have very few iterations, so do have a few (totally empirical)
                    k_n_star = (std::max)(k_n_star, 2 * m);
                }
            }

            ++iterations_;
            if (debug_verbosity_level > 1)
                SM_DEBUG ("[sm::ProgressiveSampleConsensus::computeModel] Trial %d out of %d: %d inliers (best is: %d so far).\n", iterations_, k_n_star, I_N, I_N_best);
            if (iterations_ > max_iterations_)
            {
                if (debug_verbosity_level > 0)
                    SM_DEBUG ("[sm::ProgressiveSampleConsensus::computeModel] RANSAC reached the maximum number of trials.\n");
                break;
            }
        }

        if (debug_verbosity_level > 0)
            SM_DEBUG ("[sm::ProgressiveSampleConsensus::computeModel] Model: %zu size, %d inliers.\n", model_.size (), I_N_best);

        if (model_.empty ())
        {
            inliers_.clear ();
            return (false);
        }

        // Get the set of inliers that correspond to the best model found so far
        //sac_model_->selectWithinDistance (model_coefficients_, threshold_, inliers_);
        return (true);
    }

    
} // namespace aslam
