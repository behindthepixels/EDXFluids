
#include "Fluid.h"
#include "Math/EDXMath.h"
#include "Math/Vector.h"

#include "RNG/Random.h"

#include <ppl.h>
using namespace Concurrency;

#include <list>

namespace EDX
{
	namespace FluidSim
	{
		template<uint Dim>
		float Lerp(const float val[1], const Vec<Dim, float>& vLin)
		{
			return 0.0f;
		}
		template<>
		float Lerp<1>(const float val[2], const Vec<1, float>& vLin)
		{
			return Math::Lerp(val[0], val[1], vLin[0]);
		}
		template<>
		float Lerp<2>(const float val[4], const Vec<2, float>& vLin)
		{
			return Math::BiLerp(val[0], val[1], val[2], val[3], vLin.x, vLin.y);
		}
		template<>
		float Lerp<3>(const float val[8], const Vec<3, float>& vLin)
		{
			return Math::TriLerp(val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7], vLin.x, vLin.y, vLin.z);
		}

		template<uint Dim>
		float Lerp2ndOrder(const float val[1], const Vec<Dim, float>& vLin)
		{
			return 0.0f;
		}
		template<>
		float Lerp2ndOrder<1>(const float val[4], const Vec<1, float>& vLin)
		{
			return Math::MonoCubicLerp(val[0], val[1], val[2], val[3], vLin[0]);
		}
		template<>
		float Lerp2ndOrder<2>(const float val[16], const Vec<2, float>& vLin)
		{
			return Math::MonoCubicLerp(
				Math::MonoCubicLerp(val[0], val[1], val[2], val[3], vLin[0]),
				Math::MonoCubicLerp(val[4], val[5], val[6], val[7], vLin[0]),
				Math::MonoCubicLerp(val[8], val[9], val[10], val[11], vLin[0]),
				Math::MonoCubicLerp(val[12], val[13], val[14], val[15], vLin[0]),
				vLin[1]);
		}
		template<>
		float Lerp2ndOrder<3>(const float val[64], const Vec<3, float>& vLin)
		{
			return Math::MonoCubicLerp(
				Lerp2ndOrder<2>(&val[0], Vector2(vLin.x, vLin.y)),
				Lerp2ndOrder<2>(&val[16], Vector2(vLin.x, vLin.y)),
				Lerp2ndOrder<2>(&val[32], Vector2(vLin.x, vLin.y)),
				Lerp2ndOrder<2>(&val[48], Vector2(vLin.x, vLin.y)),
				vLin[2]);
		}

		template<uint Dim>
		Vector3 Curl(const Vector3 vDvdp[1])
		{
			return 0.0f;
		}
		template<>
		Vector3 Curl<2>(const Vector3 vDvdp[2])
		{
			return Vector3(0.0f, 0.0f, Math::Curl(Vector2(vDvdp[0].x, vDvdp[0].y), Vector2(vDvdp[1].x, vDvdp[1].y)));
		}
		template<>
		Vector3 Curl<3>(const Vector3 vDvdp[3])
		{
			return Math::Curl(vDvdp[0], vDvdp[1], vDvdp[2]);
		}


		template<uint Dimension>
		void FluidSolver<Dimension>::Initialize(const Vec<Dimension, uint>& fluidDim, const Vec<Dimension, float> vFluidExtent)
		{
			// Set up simulation domain
			mvDim = fluidDim;
			mfDx = vFluidExtent[1] / mvDim[1];
			mfInvDx = 1.0f / mfDx;
			miFrameIdx = 0;
			mbAllowPartiallyFilled = false;

			for(uint i = 0; i < Math::Pow2<Dimension>::Value; i++)
			{
				for(uint d = 0; d < Dimension; d++)
					mavOffsetTable[i][d] = (i & (1 << d)) != 0;
			}
			for(uint i = 0; i < Math::Pow4<Dimension>::Value; i++)
			{
				for(uint d = 0; d < Dimension; d++)
					mavOffsetTable2ndOrder[i][d] = ((i & (3 << 2 * d)) >> 2 * d) - 1;
			}

			// Set up mac-grid
			for(auto i = 0; i < Dimension; i++)
			{
				mvExtent[i] = vFluidExtent[i];

				const Vec<Dimension, uint> velSize = mvDim + Vec<Dimension, uint>::UNIT[i];
				mVelocity[i].Init(velSize);
				mTempVel[i].Init(velSize);
				mVelocityMCBuffer[i].Init(velSize);
				mDeltaVel[i].Init(velSize);
				mVolume[i].Init(velSize);
				mFluidVolume[i].Init(velSize);

				// Init average velocity at cell center
				mAvgVelocity[i].Init(mvDim);
			}

			
			mLevelSet.Initialize(fluidDim);
			mSolidPhi.Initialize(fluidDim);
			mMarkers.Init(fluidDim);
			mPressure.Init(fluidDim);
			mDivergence.Init(fluidDim);

			mMatrix.resize(mPressure.LinearSize());

			// Initialize MacChormack buffer
			mMCBuffer.Init(fluidDim);

			// Vorticity Confinement buffer
			mCurl.Init(fluidDim);
			mCurlMag.Init(fluidDim);
			mVCForce.Init(fluidDim);

			mVisual.Init(fluidDim);
			{
				// Set all interior cells as air to calculate solid phi first
				SetMarkers([this](const Vec<Dimension, float>& vPos) -> CellType
				{
					return Air;
				});
				mSolidPhi.ReInitSDF(mMarkers, Solid);
				ComputeVolume(); 

				SetMarkers([this](const Vec<Dimension, float>& vPos) -> CellType
				{
					Vec<Dimension, float> vCenter = Vec<Dimension, float>(mvDim) * 0.5f;
					vCenter[1] = 0.55f * mvDim[1];
					if(Math::Distance(vPos, vCenter) < 0.10f * mvDim[1] || vPos.y < mvDim[1] * 0.22f)
						return Fluid;
					else
						return Air;
				});
				SetBoundaryMarkers();
			}
			//{
 		//		// Set all interior cells as air to calculate solid phi first
 		//		SetMarkers([this](const Vec<Dimension, float>& vPos) -> CellType
 		//		{
 		//			Vec<Dimension, float> vCenter = Vec<Dimension, float>(mvDim) * 0.5f;
			//		vCenter[1] = 0.5f * mvDim[1];
			//		//vCenter[2] = vPos[2];
 		//			if(Math::Distance(vPos, vCenter) < 0.2f * mvDim[1])
 		//				return Solid;
 		//			//else if(vPos[1] == 1.0f || vPos[1] == mvDim[1] - 2)
 		//			//	return Air;
 		//			else
 		//				return Fluid;
 		//		});
 		//		mSolidPhi.InitSDF([this](const Vec<Dimension, float>& vPos) -> float
 		//		{
 		//			Vec<Dimension, float> vCenter = Vec<Dimension, float>(mvDim) * 0.5f;
			//		vCenter[1] = 0.5f * mvDim[1];
			//		//vCenter[2] = vPos[2];
 		//			return Math::Distance(vPos, vCenter) - 0.2f * mvDim[1];
 		//		});
 		//		mLevelSet.InitSDF([this](const Vec<Dimension, float>& vPos) -> float
 		//		{
 		//			Vec<Dimension, float> vCenter = Vec<Dimension, float>(mvDim) * 0.5f;
			//		vCenter[1] = 0.5f * mvDim[1];
			//		//vCenter[2] = vPos[2];
 		//			return Math::Max(-(Math::Distance(vPos, vCenter) - 0.2f * mvDim[1]), Math::Max(1.0f - vPos[1], vPos[1] - (mvDim[1] - 2)));
 		//		});
 
 		//		//mSolidPhi.ReInitSDF(mMarkers, Solid);
 		//		//mLevelSet.ReInitSDF(mMarkers, Fluid);
 		//		ComputeVolume(); 
 		//		SetBoundaryMarkers();
 
 		//		for(auto dim = 0; dim < Dimension; dim++)
 		//		{
 		//			parallel_for(0, (int)mFluidVolume[dim].LinearSize(), [&](int i)
 		//			{
 		//				Vec<Dimension, uint> vIdx = mFluidVolume[dim].Index(i);
 		//				if(vIdx[dim] <= 1 || vIdx[dim] >= mLevelSet.GetPhi().Size(dim) - 1)
 		//				{
 		//					mFluidVolume[dim][vIdx] = 0.0f;
 		//					return;
 		//				}
 		//				for(auto d = 0; d < Dimension; d++)
 		//				{
 		//					if(d == dim)
 		//						continue;
 		//					if(vIdx[d] <= 0 || vIdx[d] >= mLevelSet.GetPhi().Size(d) - 1)
 		//					{
 		//						mFluidVolume[dim][vIdx] = 0.0f;
 		//						return;
 		//					}
 		//				}
 
 		//				const Vec<Dimension, uint> vOffset = Vec<Dimension, uint>::UNIT[dim];
 		//				mFluidVolume[dim][i] = Math::Clamp(FractionInside(mLevelSet.GetPhi()[vIdx - vOffset], mLevelSet.GetPhi()[vIdx]), 0.0f, 1.0f);
 		//			});
 		//		}

			//	mbRebuildPossionPerFrame = false;
 		//	}

			// Form poisson
			if(!mbRebuildPossionPerFrame)
				FormMatrixAndPreconditioner();

			if(mbUseFLIP)
				SeedParticles();

			mParticleCounts.Init(fluidDim);
			mWeights.Init(mvDim + Vec<Dimension, uint>::UNIT_SCALE);
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::SeedParticles()
		{
			// Init particles
			mParticles.clear();

			RandomGen random;
			for(auto i = 0; i < mMarkers.LinearSize(); i++)
			{
				if(!(mMarkers[i] & Fluid))
					continue;

				// Insert 2^Dim particles per cell
				const Vec<Dimension, float> vPosBase = Vec<Dimension, float>(mMarkers.Index(i))
					- 0.5f * Vec<Dimension, float>::UNIT_SCALE/* * mfDx*/;
				for(auto j = 0; j < Math::Pow2<Dimension>::Value; j++)
				{
					Vec<Dimension, float> vJitter;
					for(auto d = 0; d < Dimension; d++)
						vJitter[d] = random.Float();

					const Vec<Dimension, float> vPos = vPosBase + 0.5f * (Vec<Dimension, float>(mavOffsetTable[j]) + vJitter)/* * mfDx*/;
					mParticles.push_back(Particle<Dimension>(vPos, GetValueFace(vPos, mVelocity, true)));
				}
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::ReseedParticles()
		{
			RandomGen random;
			mParticleCounts.Clear();
			std::list<Particle<Dimension>> particleList(mParticles.begin(), mParticles.end());
			for(std::list<Particle<Dimension>>::const_iterator it = particleList.begin(); it != particleList.end(); ++it)
			{
				const Vec<Dimension, int> vIdx = Math::RoundToInt(it->vPos);

				if(mMarkers[vIdx] & Solid)
				{
					it = particleList.erase(it);
					continue;
				}

				//const Vec<Dimension, int> vIdx = it->vPos + 0.5f * Vec<Dimension, float>::UNIT_SCALE;
				mParticleCounts[vIdx]++;
			}

			mParticles.clear();
			std::copy(particleList.begin(), particleList.end(), std::back_inserter(mParticles));
			for(auto i = 0; i < mParticleCounts.LinearSize(); i++)
			{
				const Vec<Dimension, int> vIdx = mParticleCounts.Index(i);
				if(!(mMarkers[vIdx] & Fluid))
					continue;

				if(mParticleCounts[i] < Math::Pow2<Dimension>::Value)
				{
					const Vec<Dimension, float> vPosBase = Vec<Dimension, float>(vIdx)
						- 0.5f * Vec<Dimension, float>::UNIT_SCALE/* * mfDx*/;

					while(mParticleCounts[i] < Math::Pow2<Dimension>::Value)
					{
						Vec<Dimension, float> vPos = vPosBase;
						for(auto d = 0; d < Dimension; d++)
							vPos[d] += random.Float();

						mParticles.push_back(Particle<Dimension>(vPos, GetValueFace(vPos, mVelocity, true)));
						mParticleCounts[i]++;
					}
				}
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::AddSources(const float fDt)
		{
			// Add sources
			for(auto i = 0; i < mVelocity[1].LinearSize(); i++)
			{
				auto vIdx = mVelocity[1].Index(i);
				if(vIdx[1] == 2 || vIdx[1] == mvDim[1] - 2)
					mVelocity[1][i] = 0.15f * mfInvDx;
			}

		}

		template<uint Dimension>
		void FluidSolver<Dimension>::Advance(const float fDt)
		{
 			float fT = 0.0f;
 			float fStepT;
 			bool bFinished = false;
 			while(!bFinished)
 			{
 				fStepT = GetCFL();
 				if(mbUseFLIP)
 					fStepT *= 2.0f;
 
 				if(fT + fStepT >= fDt)
 				{
 					fStepT = fDt - fT;
 					bFinished = true;
 				}
 				else if(fT + 1.5f * fStepT >= fDt)
 				{
 					fStepT = 0.5f * (fDt - fT);
 				}
 
 				Step(fStepT);
 				fT += fStepT;
 			}
			
			miFrameIdx++;
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::Step(const float fDt)
		{
			if(mbUseFLIP)
			{
				if(mbReseedParticles)
					ReseedParticles();

				// For FLIP/PIC
				GetVelocityUpdate();
				TransferGridToParticles();
				AdvectParticles(fDt);
				TransferParticlesToGrid();

				// For FLIP
				SaveVelocity();
			}
			else // Semi-Lag advection with 2nd order interpolation
			{
				AdvectVelocityMacCormack(fDt, mVelocity, mTempVel);
				//for(auto dim = 0; dim < Dimension; dim++)
				//	DiffuseValueGaussSedel(fDt, 0.00038f, mVelocity[dim], mTempVel[dim]);
			}

			if(mbUseVorticityConfinement)
				VorticityConfinement(fDt, 0.3f);

			AddSources(fDt);

			UpdateMarkers();

			// Incompressibility
			Project(fDt);
			ConstrainVelocity();
			
			CalcVelocityAtCellCenter();
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::AdvectVelocitySemiLag(const float fDt,
			Array<Dimension, float> grid[Dimension],
			Array<Dimension, float> prevVelGrid[Dimension],
			Array<Dimension, float> velGrid[Dimension],
			const bool b2ndOrderLerp)
		{
			for(auto iDim = 0; iDim < Dimension; iDim++)
				swap(grid[iDim], prevVelGrid[iDim]);

			for(auto iDim = 0; iDim < Dimension; iDim++)
			{
				Vec<Dimension, float> vOffset;
				vOffset[iDim] = -0.5f;
				parallel_for(0, (int)grid[iDim].LinearSize(), [&](int i)
				{
					const Vec<Dimension, float> vPos = (Vec<Dimension, float>(grid[iDim].Index(i)) + vOffset)/* * mfDx*/;
					const Vec<Dimension, float> vVel = GetValueFace(vPos, velGrid, b2ndOrderLerp);
					if(vVel == Vec<Dimension, float>::ZERO)
						return;

					const Vec<Dimension, float> vNewPos = vPos - vVel * fDt;
					const float fNewVel = GetValueFace(vNewPos, prevVelGrid, Component(iDim), b2ndOrderLerp);
					grid[iDim][i] = fNewVel;
				});
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::AdvectVelocityMacCormack(const float fDt,
			Array<Dimension, float> velGrid[Dimension],
			Array<Dimension, float> prevVelGrid[Dimension])
		{
			AdvectVelocitySemiLag(fDt, velGrid, prevVelGrid, prevVelGrid);
			AdvectVelocitySemiLag(-fDt, velGrid, mVelocityMCBuffer, prevVelGrid);

			for(auto iDim = 0; iDim < Dimension; iDim++)
			{
				Vec<Dimension, float> vOffset;
				vOffset[iDim] = -0.5f;
				parallel_for(0, (int)velGrid[iDim].LinearSize(), [&](int i)
				{
					const Vec<Dimension, uint> vIdx = velGrid[iDim].Index(i);
					Vec<Dimension, int> vMkIdx = vIdx;
					vMkIdx[iDim] = Math::Clamp(vMkIdx[iDim], 0, mMarkers.Size(iDim) - 1);
					if((mMarkers[vMkIdx] & NearBoundary) != 0)
					{
						velGrid[iDim][i] = mVelocityMCBuffer[iDim][i];
						return;
					}

					float fUnclampedVal = mVelocityMCBuffer[iDim][i] + 0.5f * (prevVelGrid[iDim][i] - velGrid[iDim][i]);

					// Clamp the MacChormack
					Vec<Dimension, float> vPos = Vec<Dimension, float>(vIdx) + vOffset/* * mfDx*/;
					const Vec<Dimension, float> vVel = GetValueFace(vPos, prevVelGrid);

					const Vec<Dimension, float> vNewPos = vPos - vVel * fDt;
					const Vec<Dimension, int> vNewIdx = Math::FloorToInt(vNewPos - vOffset/* * mfInvDx*/);

					float fMin = Math::EDX_INFINITY, fMax = Math::EDX_NEG_INFINITY;
					for(auto j = 0; j < Math::Pow2<Dimension>::Value; j++)
					{
						Vec<Dimension, uint> vNeighborIdx = vNewIdx + mavOffsetTable[j];

						for(auto d = 0; d < Dimension; d++)
							vNeighborIdx[d] = Math::Clamp(vNeighborIdx[d], 0U, velGrid[iDim].Size(d) - 1);

						float fVal = prevVelGrid[iDim][vNeighborIdx];
						fMin = Math::Min(fMin, fVal);
						fMax = Math::Max(fMax, fVal);
					}

					velGrid[iDim][i] = Math::Clamp(fUnclampedVal, fMin, fMax);
				});
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::AdvectValueSemiLag(const float fDt,
			Array<Dimension, float>& grid,
			Array<Dimension, float>& prevGrid,
			const bool b2ndOrderLerp)
		{
			swap(grid, prevGrid);

			parallel_for(0, (int)grid.LinearSize(), [&](int i)
			{
				const Vec<Dimension, float> vPos = Vec<Dimension, float>(grid.Index(i))/* * mfDx*/;
				const Vec<Dimension, float> vVel = GetValue(vPos, mAvgVelocity);
				if(vVel == Vec<Dimension, float>::ZERO)
					return;

				const Vec<Dimension, float> vNewPos = vPos - vVel * fDt;
				const float fNewVal = GetValue(vNewPos, prevGrid, b2ndOrderLerp);
				grid[i] = fNewVal;
			});
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::AdvectValueMacCormack(const float fDt,
			Array<Dimension, float>& grid,
			Array<Dimension, float>& prevGrid)
		{
			AdvectValueSemiLag(fDt, grid, prevGrid);
			AdvectValueSemiLag(-fDt, grid, mMCBuffer);

			parallel_for(0, (int)grid.LinearSize(), [&](int i)
			{
				const Vec<Dimension, uint> vIdx = grid.Index(i);
				if((mMarkers[vIdx] & NearBoundary) != 0)
				{
					grid[i] = mMCBuffer[i];
					return;
				}

				float fUnclampedVal = mMCBuffer[i] + 0.5f * (prevGrid[i] - grid[i]);

				// Clamp the MacChormack
				const Vec<Dimension, float> vPos = Vec<Dimension, float>(vIdx)/* * mfDx*/;
				const Vec<Dimension, float> vVel = GetValueFace(vPos, mVelocity);

				const Vec<Dimension, float> vNewPos = vPos - vVel * fDt;
				const Vec<Dimension, int> vNewIdx = Math::FloorToInt(vNewPos/* * mfInvDx*/);

				float fMin = Math::EDX_INFINITY, fMax = Math::EDX_NEG_INFINITY;
				for(auto j = 0; j < Math::Pow2<Dimension>::Value; j++)
				{
					Vec<Dimension, uint> vNeighborIdx = vNewIdx + mavOffsetTable[j];

					for(auto d = 0; d < Dimension; d++)
						vNeighborIdx[d] = Math::Clamp(vNeighborIdx[d], 0U, grid.Size(d) - 1);

					float fVal = prevGrid[vNeighborIdx];
					fMin = Math::Min(fMin, fVal);
					fMax = Math::Max(fMax, fVal);
				}

				grid[i] = Math::Clamp(fUnclampedVal, fMin, fMax);
			});
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::DiffuseValueGaussSedel(const float fDt,
			const float fRate,
			Array<Dimension, float>& grid,
			Array<Dimension, float>& prevGrid)
		{
			if(fRate == 0.0f)
				return;

			float fAlpha = fDt * fRate / (mfDx * mfDx);
			float fInvEps = 1.0f / (1.0f + 2.0f * Dimension * fAlpha);

			swap(grid, prevGrid);

			for(auto k = 0; k < 20; k++)
			{
				parallel_for(0, (int)grid.LinearSize(), [&](int i)
				{
					const Vec<Dimension, uint> vIdx = grid.Index(i);
					if(GetValueFace(vIdx, mFluidVolume) == 0.0f)
						return;

					float fSum = 0.0f;
					for(auto dim = 0; dim < Dimension; dim++)
					{
						Vec<Dimension, uint> vNeighborIdx = vIdx;

						vNeighborIdx[dim] = vIdx[dim] + 1;
						fSum += grid[vNeighborIdx];

						vNeighborIdx[dim] = vIdx[dim] - 1;
						fSum += grid[vNeighborIdx];
					}

					grid[i] = (prevGrid[i] + fAlpha * fSum) * fInvEps;
				});
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::VorticityConfinement(const float fDt, const float fRate)
		{
			// Calculate Curl
			parallel_for(0, (int)mCurl.LinearSize(), [&](int i)
			{
				Vec<Dimension, float> vPos = mCurl.Index(i);
				if(!(mMarkers[vPos] & Fluid))
					return;

				Vector3 vDvdp[Dimension];
				for(auto d = 0; d < Dimension; d++)
				{
					Vec<Dimension, float> vOffset;
					vOffset[d] = 0.5f;
					vDvdp[d] = Math::ToVec3(GetValueFace(vPos + vOffset, mVelocity) - GetValueFace(vPos - vOffset, mVelocity));
				}

				mCurl[i] = Curl<Dimension>(vDvdp) * mfInvDx;
				mCurlMag[i] = Math::Length(mCurl[i]);
			});

			if(Dimension == 2)
			{
				parallel_for(0, (int)mCurl.LinearSize(), [&](int i)
				{
					mVisual[i] = mCurl[i].z;
				});
			}

			parallel_for(0, (int)mVCForce.LinearSize(), [&](int i)
			{
				Vec<Dimension, float> vPos = mVCForce.Index(i);
				if(!(mMarkers[vPos] & Fluid))
					return;

				const Vector3 vCurlNormal = Math::ToVec3(InterpolateGradient(vPos, mCurlMag) * mfInvDx);
				float fMag = Math::Length(vCurlNormal);
				if(fMag < float(Math::EDX_EPSILON))
					return;

				mVCForce[i] = Math::Cross(mCurl[i], vCurlNormal / fMag) * mfDx * fRate; 
			});

			for(auto d = 0; d < Dimension; d++)
			{
				parallel_for(0, (int)mVelocity[d].LinearSize(), [&](int i)
				{
					Vec<Dimension, uint> vIdx = mVelocity[d].Index(i);
					if(mFluidVolume[d][i] == 0.0f)
						return;

					const Vec<Dimension, int> vUpIdx = vIdx - Vec<Dimension, int>::UNIT[d];
					mVelocity[d][i] -= 0.5f * (mVCForce[vUpIdx][d] + mVCForce[vIdx][d]) * fDt/* * mfDx * mfInvDx*/;
				});
			}
		}


		template<uint Dimension>
		void FluidSolver<Dimension>::CalcVelocityAtCellCenter()
		{
			parallel_for(0, (int)mAvgVelocity[0].LinearSize(), [&](int i)
			{
				for(auto d = 0; d < Dimension; d++)
				{
					Vec<Dimension, float> vPos = mAvgVelocity[d].Index(i);
					mAvgVelocity[d][i] = GetValueFace(vPos, mVelocity, Component(d), true)/* * mfDx*/;
				}
			});
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::ConstrainVelocity()
		{
			for(auto iDim = 0; iDim < Dimension; iDim++)
				swap(mVelocity[iDim], mTempVel[iDim]);

			for(auto iDim = 0; iDim < Dimension; iDim++)
			{
				Vec<Dimension, float> vOffset;
				vOffset[iDim] = -0.5f;
				parallel_for(0, (int)mVolume[iDim].LinearSize(), [&](int i)
				{
					mVelocity[iDim][i] = mTempVel[iDim][i];

					if(mVolume[iDim][i] == 0.0f && mFluidVolume[iDim][i] > 0.0f)
					{
						const Vec<Dimension, float> vIdx = mTempVel[iDim].Index(i);
						const Vec<Dimension, float> vPos = vIdx + vOffset/* * mfDx*/;
						const Vec<Dimension, float> vVel = GetValueFace(vPos, mTempVel);
						const Vec<Dimension, float> vGrad = InterpolateGradient(vIdx, mSolidPhi.GetPhi());
						if(vGrad == Vec<Dimension, float>::ZERO)
							return;

						const Vec<Dimension, float> vNormal = Math::Normalize(vGrad);
						const float fPerpComp = Math::Dot(vVel, vNormal);
						const Vec<Dimension, float> vNewVel = vVel - fPerpComp * vNormal;
						mVelocity[iDim][i] = vNewVel[iDim];
					}
				});
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::AdvectParticles(const float fDeltaT)
		{
			const uint SUB_STEP = 5;
			const float fDt = fDeltaT * 0.2f;
			const float fBound = 0.51f;

			for(auto t = 0; t < SUB_STEP; t++)
			{
				parallel_for(0, (int)mParticles.size(), [&](int i)
				{
					Particle<Dimension>& particle = mParticles[i];

					const Vec<Dimension, float> vVel = GetValueFace(particle.vPos, mVelocity);
					if(vVel == Vec<Dimension, float>::ZERO)
						return;

					Vec<Dimension, float> vMidPt = particle.vPos + vVel * 0.5f * fDt;
					for(auto dim = 0; dim < Dimension; dim++)
						vMidPt[dim] = Math::Clamp(vMidPt[dim], fBound, float(mvDim[dim] - 1) - fBound);

					const Vec<Dimension, float> vMidVel = GetValueFace(vMidPt, mVelocity);
					particle.vPos += vMidVel * fDt;
					for(auto dim = 0; dim < Dimension; dim++)
						particle.vPos[dim] = Math::Clamp(particle.vPos[dim], fBound, float(mvDim[dim] - 1) - fBound);
				});
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::TransferParticlesToGrid()
		{
			for(auto dim = 0; dim < Dimension; dim++)
			{
				mVelocity[dim].Clear();
				mWeights.Clear();
				for(const auto& particle : mParticles)
				{
					const Vec<Dimension, float> vStaggeredPos = particle.vPos + 0.5f * Vec<Dimension, float>::UNIT[dim]/* * mfInvDx*/;
					//vStaggeredPos[dim] += 0.5f;

					Vec<Dimension, int> vBaseIdx;
					Vec<Dimension, float> vLin;
					for(auto d = 0; d < Dimension; d++)
					{
						vBaseIdx[d] = Math::FloorToInt(vStaggeredPos[d]);
						vLin[d] = vStaggeredPos[d] - vBaseIdx[d];
					}

					for(uint i = 0; i < Math::Pow2<Dimension>::Value; i++)
					{
						const Vec<Dimension, int>& vOffset = mavOffsetTable[i];

						Vec<Dimension, float> vWeightVec;
						for(auto d = 0; d < Dimension; d++)
							vWeightVec[d] = vOffset[d] == 0 ? 1.0f - vLin[d] : vLin[d];
						const float fWeight = vWeightVec.Product();
						mVelocity[dim][vBaseIdx + vOffset] += fWeight * particle.vVel[dim];
						mWeights[vBaseIdx + vOffset] += fWeight;
					}
				}

				for(auto i = 0; i < mVelocity[dim].LinearSize(); i++)
				{
					if(mVelocity[dim][i] != 0.0f)
						mVelocity[dim][i] /= mWeights[mVelocity[dim].Index(i)];
				}
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::TransferGridToParticles()
		{
			parallel_for(0, (int)mParticles.size(), [&](int i)
			{
				Particle<Dimension>& particle = mParticles[i];
				const Vec<Dimension, float> vVelPIC = GetValueFace(particle.vPos, mVelocity);
				const Vec<Dimension, float> vVelFLIP = particle.vVel + GetValueFace(particle.vPos, mDeltaVel);
				particle.vVel = Math::Lerp(vVelPIC, vVelFLIP, 0.99f);
			});
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::SaveVelocity()
		{
			for(auto dim = 0; dim < Dimension; dim++)
				mDeltaVel[dim] = mVelocity[dim];
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::GetVelocityUpdate()
		{
			for(auto dim = 0; dim < Dimension; dim++)
			{
				parallel_for(0, (int)mVelocity[dim].LinearSize(), [&](int i)
				{
					mDeltaVel[dim][i] = mVelocity[dim][i] - mDeltaVel[dim][i];
				});
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::ComputeVolume()
		{
			for(auto dim = 0; dim < Dimension; dim++)
			{
				parallel_for(0, (int)mVolume[dim].LinearSize(), [&](int i)
				{
					Vec<Dimension, uint> vIdx = mVolume[dim].Index(i);
					if(vIdx[dim] <= 1 || vIdx[dim] >= mSolidPhi.GetPhi().Size(dim) - 1)
					{
						mVolume[dim][vIdx] = 0.0f;
						return;
					}
					for(auto d = 0; d < Dimension; d++)
					{
						if(d == dim)
							continue;
						if(vIdx[d] <= 0 || vIdx[d] >= mSolidPhi.GetPhi().Size(d) - 1)
						{
							mVolume[dim][vIdx] = 0.0f;
							return;
						}
					}

					const Vec<Dimension, uint> vOffset = Vec<Dimension, uint>::UNIT[dim];
					mVolume[dim][i] = Math::Clamp(1.0f - FractionInside(mSolidPhi.GetPhi()[vIdx - vOffset], mSolidPhi.GetPhi()[vIdx]), 0.0f, 1.0f);
				});
			}
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::SetBoundaryMarkers()
		{
			parallel_for(0, (int)mMarkers.LinearSize(), [&](int i)
			{
				Vec<Dimension, int> vIdx = mMarkers.Index(i);

				for(auto d = 0; d < Dimension; d++)
				{
					for(auto i = -5; i <= 5; i++)
					{
						Vec<Dimension, int> vOffsetIdx = vIdx + i * Vec<Dimension, int>::UNIT[d];
						//vOffsetIdx[d] += i;
						vOffsetIdx[d] = Math::Clamp(vOffsetIdx[d], 0, mMarkers.Size(d) - 1);
						if(mMarkers[vOffsetIdx] & (Solid | Air))
						{
							mMarkers[vIdx] = CellType(mMarkers[vIdx] | NearBoundary);
							return;
						}
					}
				}
			});

		}

		template<uint Dimension>
		void FluidSolver<Dimension>::FormMatrixAndPreconditioner()
		{
			mMatrix.zero();
			const float fConst = 1.0f / (mfDx * mfDx);

			for(auto iLinearIdx = 0; iLinearIdx < mPressure.LinearSize(); iLinearIdx++)
			{
				Vec<Dimension, uint> vIndex = mPressure.Index(iLinearIdx);
				mDivergence[iLinearIdx] = 0;
				
				if(mbAllowPartiallyFilled ? GetValueFace(vIndex, mFluidVolume) == 0.0f : !(mMarkers[iLinearIdx] & Fluid))
					continue;

				for(auto dim = 0; dim < Dimension; dim++)
				{
					if(mbAllowPartiallyFilled || !(mMarkers[iLinearIdx - mPressure.Stride(dim)] & Solid))
					{
						const float fLeftTerm = mVolume[dim][vIndex] * fConst;
						mMatrix.add_to_element(iLinearIdx, iLinearIdx, mMarkers[iLinearIdx - mPressure.Stride(dim)] & Air ? fLeftTerm / Math::Max(mFluidVolume[dim][vIndex], 0.01f) : fLeftTerm);
						if(!(mMarkers[iLinearIdx - mPressure.Stride(dim)] & Air))
							mMatrix.add_to_element(iLinearIdx, iLinearIdx - mPressure.Stride(dim), -fLeftTerm);
					}
						
					if(mbAllowPartiallyFilled || !(mMarkers[iLinearIdx + mPressure.Stride(dim)] & Solid))
					{
						Vec<Dimension, uint> vOffsetIdx = vIndex + Vec<Dimension, uint>::UNIT[dim];
						const float fRightTerm = mVolume[dim][vOffsetIdx] * fConst;
						mMatrix.add_to_element(iLinearIdx, iLinearIdx, mMarkers[iLinearIdx + mPressure.Stride(dim)] & Air ? fRightTerm / Math::Max(mFluidVolume[dim][vOffsetIdx], 0.01f) : fRightTerm);
						if(!(mMarkers[iLinearIdx + mPressure.Stride(dim)] & Air))
							mMatrix.add_to_element(iLinearIdx, iLinearIdx + mPressure.Stride(dim), -fRightTerm);
					}
				}
			}

			mPCGSolver.form_preconditioner(mMatrix);
		}

		template<uint Dimension>
		void FluidSolver<Dimension>::Project(const float fDt)
		{
			// Matrix construction
			if(mbRebuildPossionPerFrame)
 				mMatrix.zero();

			const float fConst = 1.0f / (mfDx * mfDx);
 			const float fInvDt = 1.0f / fDt;

			parallel_for(0, (int)mPressure.LinearSize(), [&](int iLinearIdx)
			{
				Vec<Dimension, uint> vIndex = mPressure.Index(iLinearIdx); 
				mDivergence[iLinearIdx] = 0;
				
				if(mbAllowPartiallyFilled ? GetValueFace(vIndex, mFluidVolume) == 0.0f : !(mMarkers[iLinearIdx] & Fluid))
					return;

				for(auto dim = 0; dim < Dimension; dim++)
				{
					Vec<Dimension, uint> vOffsetIdx = vIndex + Vec<Dimension, uint>::UNIT[dim];
					if(mbRebuildPossionPerFrame)
					{
						if(mbAllowPartiallyFilled || !(mMarkers[iLinearIdx - mPressure.Stride(dim)] & Solid))
						{
							const float fLeftTerm = mVolume[dim][vIndex] * fConst;
							mMatrix.add_to_element(iLinearIdx, iLinearIdx, mMarkers[iLinearIdx - mPressure.Stride(dim)] & Air ? fLeftTerm / Math::Max(mFluidVolume[dim][vIndex], 0.01f) : fLeftTerm);
							if(!(mMarkers[iLinearIdx - mPressure.Stride(dim)] & Air))
								mMatrix.add_to_element(iLinearIdx, iLinearIdx - mPressure.Stride(dim), -fLeftTerm);
						}
						
						if(mbAllowPartiallyFilled || !(mMarkers[iLinearIdx + mPressure.Stride(dim)] & Solid))
						{
							const float fRightTerm = mVolume[dim][vOffsetIdx] * fConst;
							mMatrix.add_to_element(iLinearIdx, iLinearIdx, mMarkers[iLinearIdx + mPressure.Stride(dim)] & Air ? fRightTerm / Math::Max(mFluidVolume[dim][vOffsetIdx], 0.01f) : fRightTerm);
							if(!(mMarkers[iLinearIdx + mPressure.Stride(dim)] & Air))
								mMatrix.add_to_element(iLinearIdx, iLinearIdx + mPressure.Stride(dim), -fRightTerm);
						}
					}

					if(mbAllowPartiallyFilled || !(mMarkers[iLinearIdx - mPressure.Stride(dim)] & Solid))
						mDivergence[iLinearIdx] += mVolume[dim][vIndex] * mVelocity[dim][vIndex] * mfInvDx * fInvDt;
					if(mbAllowPartiallyFilled || !(mMarkers[iLinearIdx + mPressure.Stride(dim)] & Solid))
						mDivergence[iLinearIdx] -= mVolume[dim][vOffsetIdx] * mVelocity[dim][vOffsetIdx] * mfInvDx * fInvDt;
				}
			});
			
			if(mbRebuildPossionPerFrame)
				mPCGSolver.form_preconditioner(mMatrix);

			// Solve the system
			double dTolerance;
			int iIter;
			mPCGSolver.set_solver_parameters(1e-5, 1000);
			bool bSuccess = mPCGSolver.solve(mMatrix, mDivergence, mPressure, dTolerance, iIter);
			assert(bSuccess);

			// Update velocity
			for(auto dim = 0; dim < Dimension; dim++)
			{
				parallel_for(0, (int)mVelocity[dim].LinearSize(), [&](int i)
				{
					const Vec<Dimension, uint> vIdx = mVelocity[dim].Index(i);
					Vec<Dimension, uint> vOffset;
					vOffset[dim] = 1;
					if(mFluidVolume[dim][i] > 0.0f && (mbAllowPartiallyFilled || mVolume[dim][i] == 1.0f))
					{
						float fTheta = ((mMarkers[vIdx] & Air) || (mMarkers[vIdx - vOffset] & Air)) ? 1.0f / Math::Max(mFluidVolume[dim][i], 0.01f) : 1.0f;
						mVelocity[dim][i] -= fDt * (mPressure[vIdx] - mPressure[vIdx - vOffset]) * mfInvDx * fTheta;
					}
					else
						mVelocity[dim][i] = 0;
				});
			}
		}
		
		template<uint Dimension>
		float FluidSolver<Dimension>::GetCFL() const
		{
			Vec<Dimension, float> vMaxVel = Vec<Dimension, float>(Math::EDX_NEG_INFINITY);

			for(auto d = 0; d < Dimension; d++)
			{
				for(auto i = 0; i < mVelocity[d].LinearSize(); i++)
				{
					vMaxVel[d] = Math::Max(vMaxVel[d], Math::Abs(mVelocity[d][i]));
				}
			}
			
			float fCFL = Math::Length(vMaxVel);
			if(fCFL < 1e-16f)
				fCFL = 1e-16f;

			return 1.0f / fCFL;
		}
		
		void FluidSolver<2>::AddVelocitySrc(const Vec<2, int>& vIdx, const Vec<2, float>& vVel)
		{
			for(auto d = 0; d < 2; d++)
			{
				mVelocity[d][vIdx] = vVel[d] * mfInvDx;
			}
		}
		void FluidSolver<3>::AddVelocitySrc(const Vec<3, int>& vIdx, const Vec<3, float>& vVel)
		{
		}

		template<uint Dimension>
		Vec<Dimension, float> FluidSolver<Dimension>::GetValueFace(const Vec<Dimension, float>& vPos, const Array<Dimension, float> velGrid[Dimension], const bool b2ndOrderLerp) const
		{
			Vec<Dimension, float> vRet;
			Vec<Dimension, float> vIdx = vPos/* * mfInvDx*/;
			for(auto i = 0; i < Dimension; i++)
			{
				const Vec<Dimension, float> vStaggeredIdx = vIdx + 0.5f * Vec<Dimension, float>::UNIT[i];
				//vStaggeredIdx[i] += 0.5f;
				vRet[i] = b2ndOrderLerp ? InterpolateValue2ndOrder(vStaggeredIdx, velGrid[i]) : InterpolateValue(vStaggeredIdx, velGrid[i]);
			}
			return vRet;
		}

		template<uint Dimension>
		float FluidSolver<Dimension>::GetValueFace(const Vec<Dimension, float>& vPos, const Array<Dimension, float> velGrid[Dimension], const Component comp, const bool b2ndOrderLerp) const
		{
			float fRet = 0.0f;
			int iDim = int(comp);
			Vec<Dimension, float> vIdx = vPos/* * mfInvDx*/;
			const Vec<Dimension, float> vStaggeredIdx = vIdx + 0.5f * Vec<Dimension, float>::UNIT[iDim];
			//vStaggeredIdx[iDim] += 0.5f;
			fRet = b2ndOrderLerp ? InterpolateValue2ndOrder(vStaggeredIdx, velGrid[iDim]) : InterpolateValue(vStaggeredIdx, velGrid[iDim]);
			return fRet;
		}

		template<uint Dimension>
		float FluidSolver<Dimension>::GetValue(const Vec<Dimension, float>& vPos, const Array<Dimension, float>& grid, const bool b2ndOrderLerp) const
		{
			Vec<Dimension, float> vIdx = vPos/* * mfInvDx*/;
			return b2ndOrderLerp ? InterpolateValue2ndOrder(vIdx, grid) : InterpolateValue(vIdx, grid);
		}

		template<uint Dimension>
		Vec<Dimension, float> FluidSolver<Dimension>::GetValue(const Vec<Dimension, float>& vPos, const Array<Dimension, float> grid[Dimension], const bool b2ndOrderLerp) const
		{
			Vec<Dimension, float> vRet;
			Vec<Dimension, float> vIdx = vPos/* * mfInvDx*/;
			for(auto i = 0; i < Dimension; i++)
				vRet[i] = b2ndOrderLerp ? InterpolateValue2ndOrder(vIdx, grid[i]) : InterpolateValue(vIdx, grid[i]);

			return vRet;
		}


		template<uint Dimension>
		float FluidSolver<Dimension>::InterpolateValue(const Vec<Dimension, float>& vIdx, const Array<Dimension, float>& grid) const
		{
			float afVal[Math::Pow2<Dimension>::Value];

			const Vec<Dimension, int> vIdxBase = Math::FloorToInt(vIdx);
			if(vIdx == vIdxBase)
				return grid[vIdxBase];

			for(uint i = 0; i < Math::Pow2<Dimension>::Value; i++)
			{
				const Vec<Dimension, int>& vOffset = mavOffsetTable[i];
				Vec<Dimension, int> vValIdx = vIdxBase + vOffset;
				for(uint d = 0; d < Dimension; d++)
					vValIdx[d] = Math::Clamp(vValIdx[d], 0, grid.Size(d) - 1);

				afVal[i] = grid[vValIdx];
			}

			return Lerp<Dimension>(afVal, vIdx - vIdxBase);
		}

		template<uint Dimension>
		float FluidSolver<Dimension>::InterpolateValue2ndOrder(const Vec<Dimension, float>& vIdx, const Array<Dimension, float>& grid) const
		{
			float afVal[Math::Pow4<Dimension>::Value];

			const Vec<Dimension, int> vIdxBase = Math::FloorToInt(vIdx);
			if(vIdx == vIdxBase)
				return grid[vIdxBase];

			for(uint i = 0; i < Math::Pow4<Dimension>::Value; i++)
			{
				const Vec<Dimension, int>& vOffset = mavOffsetTable2ndOrder[i];
				Vec<Dimension, int> vValIdx = vIdxBase + vOffset;
				for(uint d = 0; d < Dimension; d++)
					vValIdx[d] = Math::Clamp(vValIdx[d], 0, grid.Size(d) - 1);

				afVal[i] = grid[vValIdx];
			}

			return Lerp2ndOrder<Dimension>(afVal, vIdx - vIdxBase);
		}

		template<uint Dimension>
		Vec<Dimension, float> FluidSolver<Dimension>::InterpolateGradient(const Vec<Dimension, float>& vIdx, const Array<Dimension, float>& grid) const
		{
			float afVal[Math::Pow2<Dimension>::Value];

			Vec<Dimension, uint> vIdxBase;
			for(auto d = 0; d < Dimension; d++)
				vIdxBase[d] = Math::FloorToInt(vIdx[d]);

			for(auto i = 0; i < Math::Pow2<Dimension>::Value; i++)
			{
				const Vec<Dimension, int>& vOffset = mavOffsetTable[i];
				Vec<Dimension, int> vValIdx = vIdxBase + vOffset;
				for(auto d = 0; d < Dimension; d++)
					vValIdx[d] = Math::Clamp(vValIdx[d], 0, grid.Size(d) - 1);

				afVal[i] = grid[vValIdx];
			}

			float afDif[Dimension][Math::Pow2<Dimension - 1>::Value] = {0};
			for(auto d = 0; d < Dimension; d++)
			{
				auto iStride = 1 << d, iGroup = Math::Pow2<Dimension - 1>::Value >> d;
				auto iIdx = 0;
				for(auto i = 0; i < iGroup; i++)
				{
					for(auto iLow = i * iStride * 2; iLow < i * iStride * 2 + iStride; iLow++)
						afDif[d][iIdx++] = afVal[iLow + iStride] - afVal[iLow];
				}
			}

			const Vec<Dimension, float> vLin = vIdx - vIdxBase;
			Vec<Dimension, float> vRet;
			for(auto d = 0; d < Dimension; d++)
			{
				Vec<Dimension - 1, float> vLin2;
				for(auto i = 0, idx = 0; i < Dimension; i++)
				{
					if(i != d)
						vLin2[idx++] = vLin[i];
				}
				vRet[d] = vLin != Vec<Dimension, float>::ZERO ? Lerp<Dimension - 1>(afDif[d], vLin2) : afDif[d][0];
			}

			return vRet;
		}

		template<uint Dimension>
		float FluidSolver<Dimension>::FractionInside(const float fPhiLeft, const float fPhiRight) const
		{
			if(fPhiLeft < 0.0f && fPhiRight < 0.0f)
				return 1.0f;
			else if(fPhiLeft < 0.0f && fPhiRight >= 0.0f)
				return fPhiLeft / (fPhiLeft - fPhiRight);
			else if(fPhiLeft >= 0.0f && fPhiRight < 0.0f)
				return fPhiRight / (fPhiRight - fPhiLeft);
			else
				return 0.0f;
		};

		template<uint Dimension>
		template<typename PhiFunc>
		void FluidSolver<Dimension>::SetMarkers(PhiFunc func)
		{
			for(auto i = 0; i < mMarkers.LinearSize(); i++)
			{
				Vec<Dimension, float> vPos = mMarkers.Index(i);
				mMarkers[i] = func(vPos);

				for(auto d = 0; d < Dimension; d++)
				{
					if(vPos[d] <= 0.0f || vPos[d] >= mvDim[d] - 1.0f)
					{
						mMarkers[i] = Solid;
					}
				}
			}
		}

		template class FluidSolver<2>;
		template class FluidSolver<3>;
	}
}