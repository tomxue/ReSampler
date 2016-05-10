// class for storing a single channel of waveform data in memory:

template<class T>
class channelInMemory{
public:

	// constructor:
	channelInMemory(size_t nSamples) : length(nSamples)
	{
		std::cout << nSamples << " " << length << "constructor" << std::endl;
		pData = new T(nSamples);
	}
	
	// copy constructor:
	channelInMemory(const channelInMemory<T>& C)
	{
		std::cout << "copy constructor" << std::endl;
		pData = new T(C.length);
		length = C.length;
	}

	//// Move constructor:
	//channelInMemory(channelInMemory<T>&& C)
	//{

	//}

	// Copy assignment:
	channelInMemory<T>& operator= (const channelInMemory<T>& C)
	{
		std::cout << "copy assignment" << std::endl;
		pData = new T(C.length);
		length = C.length;
	}

	void reAllocate(size_t nSamples) // this probably won't be used
	{	
		delete[] pData;
		pData = new T(nSamples);
		length = nSamples;
	}

	T& operator[](int index) 
	{
		// to-do: bounds checking
		return pData[index];
	}

	T& at(int index)
	{
		return pData[index];
	}

	~channelInMemory() 
	{
		delete[] pData;
	}
	T* pData;
private:
	
	size_t length;
};

// RAII class for storing an entire sound file in memory.
// All memory management is handled internally
//
// individual samples accessed like this: 
// soundfileInMemory snd;
// T sample = snd[channel][index]; // read
// snd[channel][index] = sample; // write

template<class T>
class soundfileInMemory{
public:
	soundfileInMemory()
	{
	}
	
	void deinterleaveFromFile(SndfileHandle* pInfile)	// read (load) file from disk and de-interleave into memory
	{
		if (pInfile->error() == SF_ERR_NO_ERROR) {
			
			// important - difference in terminology:
			// what I call "samples" , libsndfile calls "frames" (1 frame contains a sample for each channel ...)
			sf_count_t totalSamples = pInfile->frames();
			soundfileInMemory::reAllocate(pInfile->channels(), totalSamples);

			T inBuffer[BUFFERSIZE];
			size_t BufferSize = (BUFFERSIZE / nChannels) * nChannels; // round down to integer multiple of nChannels (file may have odd number of channels!)
			int sampleIndex = 0;
			sf_count_t count;
			pInfile->seek(0i64, SEEK_SET);

			do {	
				count = pInfile->read(inBuffer, BufferSize);
				for (int s = 0; s < count; s += nChannels) {
					for (int c = 0; c < nChannels; ++c) {

						channels[c].pData[s] = inBuffer[s + c];
						//	[sampleIndex] = inBuffer[s+c];
					}
					assert(sampleIndex < nSamples);
					sampleIndex++;
				}
			} while (count > 0);

		}

		else { // problem with file
			// to-do:
		}

	}
	void interleaveToFile(SndfileHandle* pOutfile);		// interleave from memory and write (save) to disk

	channelInMemory<T>& operator[](int index) 
	{
		return channels[index];
	}

	std::vector<channelInMemory<T>> channels;			// the only resource owned by this class: a vector of channels

	int nChannels;
	int nSamples;

	// To-do: keep metadata from header chunks etc

	~soundfileInMemory()
	{
		channels.clear();
	}

private:
	void reAllocate(int nCh, size_t nSamp)
	{
		channels.clear();
	//	for (int i = 0; i < nCh; ++i) {
		//	channels.insert(channels.end(), nCh, channelInMemory<T>(nSamp));
		channelInMemory<T> X(nSamp);
		channels.push_back(X);

		channelInMemory<T> Y(nSamp);
		channels.push_back(Y);

	//	}
		soundfileInMemory::nChannels = nCh;
		soundfileInMemory::nSamples = nSamp;
	}
};