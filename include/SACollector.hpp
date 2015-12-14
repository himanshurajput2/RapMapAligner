#ifndef SA_COLLECTOR_HPP
#define SA_COLLECTOR_HPP

#include "RapMapUtils.hpp"
#include "RapMapSAIndex.hpp"
#include "SASearcher.hpp"

#include <algorithm>
#include <iterator>

class SACollector {
    public:

	int SGetMax(int i, int j, char *Seq1, char *Seq2,int M , int N, int** arr, int** type,int GapPenality) {
            int Sim;
            int Similar = 10;
            int NonSimilar =3;
            int Gap = GapPenality;
            if ((char)Seq1[i-1] == (char)Seq2[j-1])
                Sim = Similar;
            else
                Sim = NonSimilar;

            int M1, M2, M3;
            M1 = arr[i - 1][j - 1] + Sim;
            M2 = arr[i][j - 1] + Gap;
            M3 = arr[i - 1][j] + Gap;
            int max = M1 >= M2 ? M1 : M2;
            int Mmax = M3 >= max ? M3 : max;

            if (Mmax == M1 && Sim == Similar)
                type[i][j] = 1;//Match
            else if (Mmax == M1 && Sim == NonSimilar)
                type[i][j] = 2;//Mismatch
            else if (Mmax == M2)
                type[i][j] = 3;//Gap in row string
            else if (Mmax == M3)
                type[i][j] = 4;//== delete for row ,Gap in column string , this is equal to a deletion of a character from row string
            return Mmax;
    }

	
    char *SAligner(char *str1, char *str2, int start1, int start2, int M, int N) {
            int **arr = NULL;
            int **type = NULL;
            char *str = NULL;
            char *cigar = NULL;
            int len = 0;
            int x=0;
            int Gap = 1;
            int count =0;
            int i =0, j =0;
            char buff[6];
			
            cigar  = (char*)(malloc((M+N+1)*sizeof(char)));

            memset(cigar,'\0',M+N+1);
    		if(0==M && 0 != N) // 1st string blank
	    	{
		        sprintf(cigar, "%d", N); 
		        strcat(cigar,"I");
		        return cigar;
		    } else if (0 == N && 0!=M) {
		        sprintf(cigar, "%d", M); 
		        strcat(cigar, "D");
			    return cigar;
		    }

            arr = (int**)(malloc((M+1)*sizeof(int*)));
            type = (int**)(malloc((M+1)*sizeof(int*)));
            str = (char*)(malloc((M+N+1)*sizeof(char)));

            memset(str,'\0',M+N+1);

            for(i =0 ;i < M+1;i++) {
                arr[i] = (int*)(malloc((N+1)*sizeof(int)));
                type[i] = (int*)(malloc((N+1)*sizeof(int)));
            }
            for(i =0 ;i < M+1;i++)
                for(int j=0; j<N+1 ;j++)
                {
                    arr[i][j] = 0;
                    type[i][j] = 0;
                }

            for(i = 0; i < N+1 ; i++)
                arr[0][i] =  i*Gap;

            for(i = 0; i < M+1 ; i++)
                arr[i][0] =  i*Gap;

            for(i = 1; i < M+1; i++)
                for ( j = 1; j < N+1; j++)
                    arr[i][j] = SGetMax(i, j, str1+start1, str2+ start2,M,N, arr, type,Gap); //gap =3

//Initialization done
            for(i = M; i>0;)
            for(int j= N ;i>0 && j>0 ;)
            {
                if (type[i][j] == 1) {
                    i--;j--;//Match
                    strcat(str,"M");
                } else if (type[i][j] == 2) {
                    i--;j--;//Mismatch
                    strcat(str,"X");
                } else if(type[i][j] == 3) {
                    j--;//gap in row
                    strcat(str,"I");
                } else if(type[i][j] == 4) {
                    i--;//Gap in column string , this is equal to a deletion of a character from row string
                    strcat(str,"D");
                }
            }

            while(i>0)
            {
                strcat(str,"D");
                i--;
            }
            while(j>0)
            {
                strcat(str,"I");
                j--;
            }

	     char *p1, *p2;

         if (!(! str || ! *str)) {
              for(p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2) {	
                    *p1 ^= *p2;
                    *p2 ^= *p1;
                    *p1 ^= *p2;
                }
	      }
        len = strlen(str);

        for( i= 0,j=0;i<len;) {
            count = 0;
            while(str[i] == str[j]) {
                count++;
                i++;
            }
	        sprintf(buff, "%d", count);
            strcat(cigar,buff);
            x = strlen(cigar);
            cigar[x++]=str[j];
            cigar[x]='\0';
            j= i;
        }
        for(int i =0 ;i < M+1;i++) {
              free(arr[i]);
              free(type[i]);
        }
              free(type);
              free(arr);
              free(str);
              return cigar;
    }


	void printSAHits(rapmap::utils::SAIntervalHit hit) {
		std::cout << "Printing Hits\n";
		std::cout << "hit.lenIn " << hit.len << std::endl;
		std::cout << "hit.queryPosIn " << hit.queryPos << std::endl;
		std::cout << std::endl;
	}
    SACollector(RapMapSAIndex* rmi) : rmi_(rmi) {}

		bool operator()(std::string& read,
				std::vector<rapmap::utils::QuasiAlignment>& hits,
				SASearcher& saSearcher,
				rapmap::utils::MateStatus mateStatus,
				bool strictCheck = false) {

        using QuasiAlignment = rapmap::utils::QuasiAlignment;
        using MateStatus = rapmap::utils::MateStatus;

        //auto& posIDs = rmi_->positionIDs;
        auto& rankDict = rmi_->rankDict;
        auto& txpStarts = rmi_->txpOffsets;
        auto& SA = rmi_->SA;
        auto& khash = rmi_->khash;
        auto& text = rmi_->seq;
        uint32_t sampFactor{1};
        auto salen = SA.size();

        auto readLen = read.length();
        auto maxDist = 1.5 * readLen;
        auto k = rapmap::utils::my_mer::k();
        auto readStartIt = read.begin();
        auto readEndIt = read.end();

        auto readRevStartIt = read.rbegin();
        auto readRevEndIt = read.rend();

        auto rb = read.begin();
        auto re = rb + k;
        int lbLeftFwd = 0, ubLeftFwd = 0;
        int lbLeftRC = 0, ubLeftRC = 0;
        int lbRightFwd = 0, ubRightFwd = 0;
        int lbRightRC = 0, ubRightRC = 0;
        int matchedLen;

        uint32_t fwdHit{0};
        uint32_t rcHit{0};

        bool foundHit = false;
        bool isRev = false;
        rapmap::utils::my_mer mer;
        rapmap::utils::my_mer rcMer;

        enum HitStatus { ABSENT = -1, UNTESTED = 0, PRESENT = 1 };
        // Record if k-mers are hits in the
        // fwd direction, rc direction or both
        struct KmerDirScore {
            KmerDirScore(rapmap::utils::my_mer kmerIn, HitStatus fwdScoreIn, HitStatus rcScoreIn) :
                kmer(kmerIn), fwdScore(fwdScoreIn), rcScore(rcScoreIn) {}
            KmerDirScore() : fwdScore(UNTESTED), rcScore(UNTESTED) {}
            rapmap::utils::my_mer kmer;
            HitStatus fwdScore{UNTESTED};
            HitStatus rcScore{UNTESTED};
        };

        // This allows implementing our heurisic for comparing
        // forward and reverse-complement strand matches
        std::vector<KmerDirScore> kmerScores;

        using rapmap::utils::SAIntervalHit;

        std::vector<SAIntervalHit> fwdSAInts;
        std::vector<SAIntervalHit> rcSAInts;

        std::vector<uint32_t> leftTxps, leftTxpsRC;
        std::vector<uint32_t> rightTxps, rightTxpsRC;
        int maxInterval{1000};

        // The number of bases that a new query position (to which
        // we skipped) should overlap the previous extension. A
        // value of 0 means no overlap (the new search begins at the next
        // base) while a value of (k - 1) means that k-1 bases (one less than
        // the k-mer size) must overlap.
        int skipOverlap = k-1;
        // Number of nucleotides to skip when encountering a homopolymer k-mer.
        int homoPolymerSkip = k/2;

        // Find a hit within the read
        // While we haven't fallen off the end
        while (re < read.end()) {

            // Get the k-mer at the current start position.
            // And make sure that it's valid (contains no Ns).
            auto pos = std::distance(readStartIt, rb);
            auto invalidPos = read.find_first_of("nN", pos);
            if (invalidPos < pos + k) {
                rb = read.begin() + invalidPos + 1;
                re = rb + k;
                continue;
            }

            // If the next k-bases are valid, get the k-mer and
            // reverse complement k-mer
            mer = rapmap::utils::my_mer(read.c_str() + pos);
            if (mer.is_homopolymer()) { rb += homoPolymerSkip; re += homoPolymerSkip; continue; }
            rcMer = mer.get_reverse_complement();

            // See if we can find this k-mer in the hash
            auto merIt = khash.find(mer.get_bits(0, 2*k));
            auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));

            // If we can find the k-mer in the hash, get its SA interval
            if (merIt != khash.end()) {
                int lb = merIt->second.begin;
                int ub = merIt->second.end;

                // lb must be 1 *less* then the current lb
                auto lbRestart = std::max(static_cast<int>(0), lb-1);
                // Extend the SA interval using the read sequence as far as
                // possible
                std::tie(lbLeftFwd, ubLeftFwd, matchedLen) =
                    saSearcher.extendSearchNaive(lbRestart, ub, k, rb, readEndIt);

                // If the SA interval is valid, and not too wide, then record
                // the hit.
                int diff = ubLeftFwd - lbLeftFwd;
                if (ubLeftFwd > lbLeftFwd and diff < maxInterval) {
                    auto queryStart = std::distance(read.begin(), rb);
                    fwdSAInts.emplace_back(lbLeftFwd, ubLeftFwd, matchedLen, queryStart, false);
                    if (strictCheck) {
                        ++fwdHit;
                        // If we also match this k-mer in the rc direction
                        if (rcMerIt != khash.end()) {
                            ++rcHit;
                            kmerScores.emplace_back(mer, PRESENT, PRESENT);
                        } else { // Otherwise it doesn't match in the rc direction
                            kmerScores.emplace_back(mer, PRESENT, ABSENT);
                        }

                        // If we didn't end the match b/c we exhausted the query
                        // test the mismatching k-mer in the other strand
                        // TODO: check for 'N'?
                        if (rb + matchedLen < readEndIt){
                            auto kmerPos = std::distance(readStartIt, rb + matchedLen - skipOverlap);
                            mer = rapmap::utils::my_mer(read.c_str() + kmerPos);
                            kmerScores.emplace_back(mer , ABSENT, UNTESTED);
                        }
                    } else { // no strict check
                        ++fwdHit;
                        if (rcMerIt != khash.end()) { ++rcHit; }
                    }
                }
            }

            // See if the reverse complement k-mer is in the hash
            if (rcMerIt != khash.end()) {
                lbLeftRC = rcMerIt->second.begin;
                ubLeftRC = rcMerIt->second.end;
                if (ubLeftRC > lbLeftRC) {
                    // The original k-mer didn't match in the foward direction
                    if (!fwdHit) {
                        ++rcHit;
                        if (strictCheck) {
                            kmerScores.emplace_back(rcMer, ABSENT, PRESENT);
                        }
                    }
                }
            }

            // If we had a hit with either k-mer then we can
            // break out of this loop to look for the next informative position
            if (fwdHit + rcHit > 0) {
                foundHit = true;
                break;
            }
            ++rb; ++re;
        }

        // If we went the entire length of the read without finding a hit
        // then we can bail.
        if (!foundHit) { return false; }

        bool lastSearch{false};
        // If we had a hit on the forward strand
        if (fwdHit) {

            // The length of this match
            int32_t matchLen = fwdSAInts.front().len;

            // The iterator to where this match began
            rb = read.begin() + fwdSAInts.front().queryPos;

            // [lb, ub) is the suffix array interval for the MMP (maximum mappable prefix)
            // of the k-mer we found.  The NIP (next informative position) in the sequence
            // is the position after the LCE (longest common extension) of
            // T[SA[lb]:] and T[SA[ub-1]:]
            auto remainingLength = std::distance(rb + matchLen, readEndIt);
            int32_t lce = saSearcher.lce(
                    lbLeftFwd, ubLeftFwd-1, matchLen, remainingLength);
            auto fwdSkip = std::max(matchLen - skipOverlap, lce - skipOverlap);

            size_t nextInformativePosition = std::min(
                    std::max(0, static_cast<int>(readLen)- static_cast<int>(k)),
                    static_cast<int>(std::distance(readStartIt, rb) + fwdSkip)
                    );

            rb = read.begin() + nextInformativePosition;
            re = rb + k;

            size_t invalidPos{0};
            while (re <= readEndIt) {
                // The offset into the string
                auto pos = std::distance(readStartIt, rb);

                // The position of the first N in the k-mer (if there is one)
                // If we have already verified there are no Ns in the remainder
                // of the string (invalidPos is std::string::npos) then we can
                // skip this test.
                if (invalidPos != std::string::npos) {
                    invalidPos = read.find_first_of("nN", pos);
                }

                // If the first N is within k bases, then this k-mer is invalid
                if (invalidPos < pos + k) {
                    // A valid k-mer can't start until after the 'N'
                    nextInformativePosition = invalidPos + 1;
                    rb = read.begin() + nextInformativePosition;
                    re = rb + k;
                    // Go to the next iteration of the while loop
                    continue;
                }

                // If the current end position is valid
                if (re <= readEndIt) {

                    mer = rapmap::utils::my_mer(read.c_str() + pos);
                    if (mer.is_homopolymer()) { rb += homoPolymerSkip; re = rb + k; continue; }
                    auto merIt = khash.find(mer.get_bits(0, 2*k));

                    if (merIt != khash.end()) {
                        if (strictCheck) {
                            ++fwdHit;
                            kmerScores.emplace_back(mer, PRESENT, UNTESTED);
                            auto rcMer = mer.get_reverse_complement();
                            auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));
                            if (rcMerIt != khash.end()) {
                                ++rcHit;
                                kmerScores.back().rcScore = PRESENT;
                            }
                        }

                      lbRightFwd = merIt->second.begin;
                      ubRightFwd = merIt->second.end;

                        // lb must be 1 *less* then the current lb
                        lbRightFwd = std::max(0, lbRightFwd - 1);
                        std::tie(lbRightFwd, ubRightFwd, matchedLen) =
                            saSearcher.extendSearchNaive(lbRightFwd, ubRightFwd,
                                    k, rb, readEndIt);

                        int diff = ubRightFwd - lbRightFwd;
                        if (ubRightFwd > lbRightFwd and diff < maxInterval) {
                            auto queryStart = std::distance(read.begin(), rb);
                            fwdSAInts.emplace_back(lbRightFwd, ubRightFwd, matchedLen, queryStart, false);
                            // If we didn't end the match b/c we exhausted the query
                            // test the mismatching k-mer in the other strand
                            // TODO: check for 'N'?
                            if (strictCheck and rb + matchedLen < readEndIt){
                                auto kmerPos = std::distance(readStartIt, rb + matchedLen - skipOverlap);
                                mer = rapmap::utils::my_mer(read.c_str() + kmerPos);
                                kmerScores.emplace_back(mer , ABSENT, UNTESTED);
                            }
                        }

                        if (lastSearch) { break; }
                        auto mismatchIt = rb + matchedLen;
                        if (mismatchIt < readEndIt) {
                            auto remainingDistance = std::distance(mismatchIt, readEndIt);
                            auto lce = saSearcher.lce(lbRightFwd, ubRightFwd-1, matchedLen, remainingDistance);

                            // Where we would jump if we just used the MMP
                            auto skipMatch = mismatchIt - skipOverlap;
                            // Where we would jump if we used the LCE
                            auto skipLCE = rb + lce - skipOverlap;
                            // Pick the larger of the two
                            rb = std::max(skipLCE, skipMatch);
                            if (rb > (readEndIt - k)) {
                                rb = readEndIt - k;
                                lastSearch = true;
                            }
                            re = rb + k;
                        } else {
                            lastSearch = true;
                            rb = readEndIt - k;
                            re = rb + k;
                        }

                    } else {
                        rb += sampFactor;
                        re = rb + k;
                    }
                }
            }
        }

        lastSearch = false;
        if (rcHit >= fwdHit) {
            size_t pos{read.length() - k};

            auto revReadStartIt = read.rend();
            auto revReadEndIt = read.rend();

            auto revRB = read.rbegin();
            auto revRE = revRB + k;

            auto invalidPosIt = revRB;
            while (revRE <= revReadEndIt){

                revRE = revRB + k;
                if (revRE > revReadEndIt) { break; }

                // See if this k-mer would contain an N
                // only check if we don't yet know that there are no remaining
                // Ns
                if (invalidPosIt != revReadEndIt) {
                    invalidPosIt = std::find_if(revRB, revRE,
                                                 [](const char c) -> bool {
                                                     return c == 'n' or c == 'N';
                                                 });
                }

                // If we found an N before the end of the k-mer
                if (invalidPosIt < revRE) {
                    // Skip to the k-mer starting at the next position
                    // (i.e. right past the N)
                    revRB = invalidPosIt + 1;
                    continue;
                }

                // The distance from the beginning of the read to the
                // start of the k-mer
                pos = std::distance(revRE, revReadEndIt);

                // Get the k-mer and query it in the hash
                mer = rapmap::utils::my_mer(read.c_str() + pos);
                if (mer.is_homopolymer()) { revRB += homoPolymerSkip; revRE += homoPolymerSkip; continue; }
                rcMer = mer.get_reverse_complement();
                auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));

                // If we found the k-mer
                if (rcMerIt != khash.end()) {
                    if (strictCheck) {
                        ++rcHit;
                        kmerScores.emplace_back(mer, UNTESTED, PRESENT);
                        auto merIt = khash.find(mer.get_bits(0, 2*k));
                        if (merIt != khash.end()) {
                            ++fwdHit;
                            kmerScores.back().fwdScore = PRESENT;
                        }
                    }


                    lbRightRC = rcMerIt->second.begin;
                    ubRightRC = rcMerIt->second.end;

                    // lb must be 1 *less* then the current lb
                    // We can't move any further in the reverse complement direction
                    lbRightRC = std::max(0, lbRightRC - 1);
                    std::tie(lbRightRC, ubRightRC, matchedLen) =
                        saSearcher.extendSearchNaive(lbRightRC, ubRightRC, k,
                                revRB, revReadEndIt, true);

                    int diff = ubRightRC - lbRightRC;
                    if (ubRightRC > lbRightRC and diff < maxInterval) {
                        auto queryStart = std::distance(read.rbegin(), revRB);
                        rcSAInts.emplace_back(lbRightRC, ubRightRC, matchedLen, queryStart, true);
                        // If we didn't end the match b/c we exhausted the query
                        // test the mismatching k-mer in the other strand
                        // TODO: check for 'N'?
                        if (strictCheck and revRB + matchedLen < revReadEndIt){
                            auto kmerPos = std::distance(revRB + matchedLen, revReadEndIt);
                            mer = rapmap::utils::my_mer(read.c_str() + kmerPos);
                            kmerScores.emplace_back(mer , UNTESTED, ABSENT);
                        }
                    }

                    if (lastSearch) { break; }
                    auto mismatchIt = revRB + matchedLen;
                    if (mismatchIt < revReadEndIt) {
                        auto remainingDistance =
                            std::distance(mismatchIt, revReadEndIt);
                        auto lce = saSearcher.lce(lbRightRC, ubRightRC-1,
                                                  matchedLen, remainingDistance);

                        // Where we would jump if we just used the MMP
                        auto skipMatch = mismatchIt - skipOverlap;
                        // Where we would jump if we used the lce
                        auto skipLCE = revRB + lce - skipOverlap;
                        // Choose the larger of the two
                        revRB = std::max(skipLCE, skipMatch);
                        if (revRB > (revReadEndIt - k)) {
                            revRB = revReadEndIt - k;
                            lastSearch = true;
                        }
                        revRE = revRB + k;

                    } else {
                        lastSearch = true;
                        revRB = revReadEndIt - k;
                        revRE = revRB + k;
                    }

                } else {
                    revRB += sampFactor;
                    revRE = revRB + k;
                }
            }
        }


        if (strictCheck) {
            // The first two conditions shouldn't happen
            // but I'm just being paranoid here
            if (fwdHit > 0 and rcHit == 0) {
                rcSAInts.clear();
            } else if (rcHit > 0 and fwdHit == 0) {
                fwdSAInts.clear();
            } else {
                // Compute the score for the k-mers we need to
                // test in both the forward and rc directions.
                int32_t fwdScore{0};
                int32_t rcScore{0};
                // For every kmer score structure
                for (auto& kms : kmerScores) {
                    // If the forward k-mer is untested, then test it
                    if (kms.fwdScore == UNTESTED) {
                        auto merIt = khash.find(mer.get_bits(0, 2*k));
                        kms.fwdScore = (merIt != khash.end()) ? PRESENT : ABSENT;
                    }
                    // accumulate the score
                    fwdScore += kms.fwdScore;

                    // If the rc k-mer is untested, then test it
                    if (kms.rcScore == UNTESTED) {
                        rcMer = kms.kmer.get_reverse_complement();
                        auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));
                        kms.rcScore = (rcMerIt != khash.end()) ? PRESENT : ABSENT;
                    }
                    // accumulate the score
                    rcScore += kms.rcScore;
                }
                // If the forward score is strictly greater
                // then get rid of the rc hits.
                if (fwdScore > rcScore) {
                    rcSAInts.clear();
                } else if (rcScore > fwdScore) {
                    // If the rc score is strictly greater
                    // get rid of the forward hits
                    fwdSAInts.clear();
                }
            }
        }

        /* VERY SIMPLE HEURISTIC */
        /*
           if (fwdHit > rcHit) {
           rcSAInts.clear();
           } else if (rcHit > fwdHit) {
           fwdSAInts.clear();
           }
           */

        auto fwdHitsStart = hits.size();

		std::string cigar;
        // If we had > 1 forward hit
        if (fwdSAInts.size() > 1) {
			auto processedHits = rapmap::hit_manager::intersectSAHits(fwdSAInts, *rmi_);
		for(std::map<int, rapmap::utils::ProcessedSAHit>::iterator iter = processedHits.begin(); iter != processedHits.end(); ++iter) {
                char *c = NULL;
				int k =  iter->first;
				cigar = "";
				std::string transcript = rmi_->seq.substr(rmi_->txpOffsets[k], rmi_->txpLens[k]);
				rapmap::utils::ProcessedSAHit val = iter->second;
				int transcriptAlignStart = 0, transcriptAlignEnd = 0, queryAlignStart = 0, queryAlignEnd = 0;
				std::vector<rapmap::utils::SATxpQueryPos> myVec = val.tqvec;
				for(int i=0;i<myVec.size();i++) {
					// Append fwdSAInts[i].len number of M's to cigar

					//std::cout << "Genome position:" << myVec[i].pos << std::endl;
					//std::cout << "Query position:" << myVec[i].queryPos << std::endl;

					transcriptAlignEnd = myVec[i].pos;
					queryAlignEnd = myVec[i].queryPos;

					// Call Aigner

					if(queryAlignEnd < queryAlignStart || transcriptAlignEnd < transcriptAlignStart) continue;
					else {
						// Add as many Number of M's as the fwdSAInts[i]
						/*std::cout << "queryAlignStart: " << queryAlignStart << std::endl;
						std::cout << "queryAlignEnd: " << queryAlignEnd << std::endl;
						std::cout << "transcriptAlignStart: " << transcriptAlignStart << std::endl;
						std::cout << "transcriptAlignEnd: " << transcriptAlignEnd<< std::endl;*/

						c = SAligner(const_cast<char*>(read.c_str()), const_cast<char*>(transcript.c_str()), queryAlignStart, transcriptAlignStart,
								queryAlignEnd - queryAlignStart, transcriptAlignEnd - transcriptAlignStart);
						cigar = cigar + std::string (c);
						cigar = cigar + std::to_string(fwdSAInts[i].len) + "M";
						transcriptAlignStart = transcriptAlignEnd + fwdSAInts[i].len;
						queryAlignStart = queryAlignEnd + fwdSAInts[i].len;
						free(c);
					}

					// We take fwdSAInts.length and add M*fwdSAInts.length to cigar string
					// We run Aligner here and update the M/D/I/S counts
				}
				iter->second.cigar_string = cigar;
				//for(int i=0;i<myVec.size();i++)
				//	printSAHits(fwdSAInts[i]);
			}
            rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist, hits, mateStatus);
        } else if (fwdSAInts.size() == 1) { // only 1 hit!
            auto& saIntervalHit = fwdSAInts.front();
            auto initialSize = hits.size();
            for (int i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
                auto globalPos = SA[i];
                auto txpID = rmi_->transcriptAtPosition(globalPos);
                // the offset into this transcript
                auto pos = globalPos - txpStarts[txpID];
                hits.emplace_back(txpID, pos, true, readLen);
                hits.back().mateStatus = mateStatus;
            }
            // Now sort by transcript ID (then position) and eliminate
            // duplicates
            auto sortStartIt = hits.begin() + initialSize;
            auto sortEndIt = hits.end();
            std::sort(sortStartIt, sortEndIt,
                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    if (a.tid == b.tid) {
                    return a.pos < b.pos;
                    } else {
                    return a.tid < b.tid;
                    }
                    });
            auto newEnd = std::unique(hits.begin() + initialSize, hits.end(),
                    [] (const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    return a.tid == b.tid;
                    });
            hits.resize(std::distance(hits.begin(), newEnd));
        }
        auto fwdHitsEnd = hits.size();

        auto rcHitsStart = fwdHitsEnd;
        // If we had > 1 rc hit
        if (rcSAInts.size() > 1) {
            auto processedHits = rapmap::hit_manager::intersectSAHits(rcSAInts, *rmi_);
            rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist, hits, mateStatus);
        } else if (rcSAInts.size() == 1) { // only 1 hit!
            auto& saIntervalHit = rcSAInts.front();
            auto initialSize = hits.size();
            for (int i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
                auto globalPos = SA[i];
                auto txpID = rmi_->transcriptAtPosition(globalPos);
                // the offset into this transcript
                auto pos = globalPos - txpStarts[txpID];
                hits.emplace_back(txpID, pos, false, readLen);
                hits.back().mateStatus = mateStatus;
            }
            // Now sort by transcript ID (then position) and eliminate
            // duplicates
            auto sortStartIt = hits.begin() + rcHitsStart;
            auto sortEndIt = hits.end();
            std::sort(sortStartIt, sortEndIt,
                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    if (a.tid == b.tid) {
                    return a.pos < b.pos;
                    } else {
                    return a.tid < b.tid;
                    }
                    });
            auto newEnd = std::unique(sortStartIt, sortEndIt,
                    [] (const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    return a.tid == b.tid;
                    });
            hits.resize(std::distance(hits.begin(), newEnd));
        }
        auto rcHitsEnd = hits.size();

        // If we had both forward and RC hits, then merge them
        if ((fwdHitsEnd > fwdHitsStart) and (rcHitsEnd > rcHitsStart)) {
            // Merge the forward and reverse hits
            std::inplace_merge(hits.begin() + fwdHitsStart, hits.begin() + fwdHitsEnd, hits.begin() + rcHitsEnd,
                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    return a.tid < b.tid;
                    });
            // And get rid of duplicate transcript IDs
            auto newEnd = std::unique(hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd,
                    [] (const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    return a.tid == b.tid;
                    });
            hits.resize(std::distance(hits.begin(), newEnd));
        }
        // Return true if we had any valid hits and false otherwise.
        return foundHit;
        }

    private:
        RapMapSAIndex* rmi_;
	/*
    std::string revComp(std::string& str) {
        std::string rc(str.size(), 'A');
        auto outIt = rc.begin();
        for (auto it = str.rbegin(); it != str.rend(); ++it, ++outIt) {
            switch (*it) {
                case 'a':
                case 'A':
                    (*outIt) = 'T';
                    break;
                case 'c':
                case 'C':
                    (*outIt) = 'G';
                    break;
                case 'g':
                case 'G':
                    (*outIt) = 'C';
                    break;
                case 't':
                case 'T':
                    (*outIt) = 'A';
                    break;
                default:
                    (*outIt) = 'A';
            }
        }
        return rc;
    }
	*/

};

#endif // SA_COLLECTOR_HPP
