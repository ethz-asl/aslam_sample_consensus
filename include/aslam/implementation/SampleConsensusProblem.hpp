namespace aslam {

    template<typename M>
    SampleConsensusProblem<M>::SampleConsensusProblem(bool randomSeed) : 
        max_sample_checks_(10)
    {
        rng_dist_.reset(new boost::uniform_int<> (0, std::numeric_limits<int>::max ()));
        // Create a random number generator object
        if (randomSeed)
            rng_alg_.seed (static_cast<unsigned> (std::time(0)));
        else
            rng_alg_.seed (12345u);

        rng_gen_.reset (new boost::variate_generator<boost::mt19937&, boost::uniform_int<> > (rng_alg_, *rng_dist_)); 

    }
    
    template<typename M>
    SampleConsensusProblem<M>::~SampleConsensusProblem()
    {

    }
    
    template<typename M>
    bool SampleConsensusProblem<M>::isSampleGood(const std::vector<int> & sample) const
    {
        // Default implementation
        return true;
    }

    template<typename M>
    void SampleConsensusProblem<M>::drawIndexSample (std::vector<int> & sample)
      {
        size_t sample_size = sample.size();
        size_t index_size = shuffled_indices_.size();
        for (unsigned int i = 0; i < sample_size; ++i)
        {
            // The 1/(RAND_MAX+1.0) trick is when the random numbers are not uniformly distributed and for small modulo
            // elements, that does not matter (and nowadays, random number generators are good)
            //std::swap (shuffled_indices_[i], shuffled_indices_[i + (rand () % (index_size - i))]);
            std::swap (shuffled_indices_[i], shuffled_indices_[i + (rnd () % (index_size - i))]);
        }
        std::copy (shuffled_indices_.begin (), shuffled_indices_.begin () + sample_size, sample.begin ());
      }


    template<typename M>
    void SampleConsensusProblem<M>::getSamples(int &iterations, std::vector<int> &samples)
    {
        // We're assuming that indices_ have already been set in the constructor
        if (indices_->size() < (size_t)getSampleSize())
        {
          fprintf(stderr,"[sm::SampleConsensusModel::getSamples] Can not select %zu unique points out of %zu!\n",
                     (size_t) getSampleSize(), indices_->size ());
          // one of these will make it stop :)
          samples.clear();
          iterations = INT_MAX - 1;
          return;
        }

        // Get a second point which is different than the first
        samples.resize( getSampleSize() );

        for (int iter = 0; iter < max_sample_checks_; ++iter)
        {
            drawIndexSample (samples);
          
          // If it's a good sample, stop here
          if (isSampleGood (samples))
            return;
        }
        fprintf(stdout,"[sm::SampleConsensusModel::getSamples] WARNING: Could not select %d sample points in %d iterations!\n", getSampleSize (), max_sample_checks_);
        samples.clear ();

    }
    
    template<typename M>
    boost::shared_ptr <std::vector<int> > SampleConsensusProblem<M>::getIndices () const
    {
        return indices_;
    }

    /** \brief Compute all distances from the cloud data to a given model. Pure .
     * 
     * \param[in] model_coefficients the coefficients of a model that we need to compute distances to 
     * \param[out] distances the resultant estimated distances
     */

    template<typename M>
    void SampleConsensusProblem<M>::getDistancesToModel (const model_t & model_coefficients, 
                                                std::vector<double> & distances)
    {
        getSelectedDistancesToModel(model_coefficients, *indices_, distances);
    }

    template<typename M>
    void SampleConsensusProblem<M>::setUniformIndices(int N)
    {
        indices_.reset( new std::vector<int>() );
        indices_->resize(N);
        for(int i = 0; i < N; ++i)
        {
            (*indices_)[i] = i;
        }
        shuffled_indices_ = *indices_;
    }

    template<typename M>
    void SampleConsensusProblem<M>::setIndices(const std::vector<int> & indices)
    {
        indices_.reset( new std::vector<int>(indices) );
        shuffled_indices_ = *indices_;
    }


    template<typename M>
    int SampleConsensusProblem<M>::rnd ()
    {
        return ((*rng_gen_) ());
    }


    template<typename M>
    void SampleConsensusProblem<M>::selectWithinDistance (const model_t & model_coefficients, 
                                                 const double threshold,
                                                 std::vector<int> &inliers)
    {
        std::vector<double> dist;
        dist.reserve(indices_->size());
        getDistancesToModel(model_coefficients, dist);

        inliers.clear();
        inliers.reserve(indices_->size());
        for(size_t i = 0; i < dist.size(); ++i)
        {
            if(dist[i] < threshold)
            {
                inliers.push_back( (*indices_)[i] );
            }
        }
    }

    template<typename M>
    int SampleConsensusProblem<M>::countWithinDistance (const model_t & model_coefficients, 
                                               const double threshold)
    {
        std::vector<double> dist;
        dist.reserve(indices_->size());
        getDistancesToModel(model_coefficients, dist);

        int count = 0;
        for(size_t i = 0; i < dist.size(); ++i)
        {
            if(dist[i] < threshold)
            {
                ++count;
            }
        }
        
        return count;
    }



} // namespace aslam
