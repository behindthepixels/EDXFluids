
#include "LevelSet.h"
#include "Math/EDXMath.h"
#include "Math/Vector.h"

#include <ppl.h>
using namespace Concurrency;

namespace EDX
{
	namespace FluidSim
	{
		template<uint Dimension>
		void LevelSet<Dimension>::Initialize(const Vec<Dimension, int>& vDim)
		{
			mPhi.Init(vDim);
			mPhiPrev.Init(vDim);
			mClosestSurfIdx.Init(vDim);
		}

		template<uint Dimension>
		void LevelSet<Dimension>::ReInitSDF(const Array<Dimension, CellType>& markers, const CellType type, const float fNarrowBand)
		{
			parallel_for(0, (int)mPhi.LinearSize(), [&](int i)
			{
				if((markers[i] & type) != 0)
					mPhi[i] = -1.0f;
				else
					mPhi[i] = 1.0f;
			});
			
			struct MarchItem
			{
				Vec<Dimension, int> vCoord;
				Vec<Dimension, float> vClosest;
				Vec<Dimension, int> vLiquidCrd;
				float fDist;

				MarchItem(const Vec<Dimension, int>& coord,
					const Vec<Dimension, float>& closest,
					const Vec<Dimension, int> liqCrd,
					float dist)
					: vCoord(coord), vClosest(closest), vLiquidCrd(liqCrd), fDist(dist) {}

				bool operator < (const MarchItem& rhs) const
				{
					return fDist > rhs.fDist;
				}
			};

			swap(mPhi, mPhiPrev);

			vector<MarchItem> marchHeap;
			for(auto i = 0; i < mPhi.LinearSize(); i++)
			{
				Vec<Dimension, float> vClosest;
				Vec<Dimension, float> vClosestCrd;
				Vec<Dimension, float> vLiquidCrd;
				float fMinDist = Math::EDX_INFINITY;

				const Vec<Dimension, int> vIdx = mPhiPrev.Index(i);
				const float fCenter = mPhiPrev[i];
				auto EdgeDetect = [&]() -> bool
				{
					bool bEdge = false;

					for(auto dim = 0; dim < Dimension; dim++)
					{
						for(auto j = -1; j <= 1; j += 2)
						{
							const Vec<Dimension, int> vNeighborIdx = vIdx + j * Vec<Dimension, int>::UNIT[dim];

							if(vNeighborIdx[dim] < 0 || vNeighborIdx[dim] >= mPhi.Size(dim))
								continue;

							float fNeighbor = mPhiPrev[vNeighborIdx];
							if((markers[vIdx] & type) != (markers[vNeighborIdx] & type))
							{
								bEdge = true;
								float fDist = Math::LinStep(0.0f, fCenter, fNeighbor);
								if(fDist < fMinDist)
								{
									fMinDist = fDist;
									vClosestCrd = vNeighborIdx;
									vLiquidCrd = fCenter <= 0.0f ? vIdx : vNeighborIdx;
								}
							}
						}
					}

					return bEdge;
				};

				// Detect initials surface boundary
				if(EdgeDetect())
				{
					const Vec<Dimension, float> vBoundary = Math::Lerp(Vec<Dimension, float>(vIdx), vClosestCrd, fMinDist);
					MarchItem item(vIdx, vBoundary, vLiquidCrd, Math::Length(vBoundary - vIdx));
					marchHeap.push_back(item);
				}

				mPhi[i] = float(Math::EDX_INFINITY);
				mClosestSurfIdx[i] = Vec<Dimension, int>::ZERO;
			}

			// Fast marching
			std::make_heap(marchHeap.begin(), marchHeap.end());
			while(!marchHeap.empty())
			{
				std::pop_heap(marchHeap.begin(), marchHeap.end());
				MarchItem item = marchHeap.back();
				if(item.fDist > fNarrowBand)
					break;
				marchHeap.pop_back();

				if(mPhi[item.vCoord] == float(Math::EDX_INFINITY))
				{
					mPhi[item.vCoord] = item.fDist;
					mClosestSurfIdx[item.vCoord] = item.vLiquidCrd;

					// Check neighbors
					for(auto dim = 0; dim < Dimension; dim++)
					{
						for(auto j = -1; j <= 1; j += 2)
						{
							const Vec<Dimension, int> vNeighborIdx = item.vCoord + j * Vec<Dimension, int>::UNIT[dim];

							if(vNeighborIdx[dim] < 0 || vNeighborIdx[dim] >= mPhi.Size(dim))
								continue;

							float fDist = Math::Distance(vNeighborIdx, item.vClosest);
							MarchItem newItem = MarchItem(vNeighborIdx, item.vClosest, item.vLiquidCrd, fDist);
							marchHeap.push_back(newItem);
							std::push_heap(marchHeap.begin(), marchHeap.end());
						}
					}
				}
			}

			for(auto i = 0; i < mPhi.LinearSize(); i++)
			{
				mPhi[i] *= !(markers[i] & type) ? 1.0f : -1.0f;
			}
		}
		
		template class LevelSet<2>;
		template class LevelSet<3>;
	}
}