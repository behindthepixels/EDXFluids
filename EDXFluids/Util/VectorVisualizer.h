#pragma once

#include "EDXPrerequisites.h"
#include "../Base/Fluid.h"
#include "Math/BoundingBox.h"
#include "Core/Random.h"


class glShader;

namespace EDX
{
	namespace FluidSim
	{
		// 2D Visualizer
		struct StreamParticle2D
		{
			struct Pos
			{
				float x;
				float y;

				Pos(float fX = float(Math::EDX_NEG_INFINITY), float fY = float(Math::EDX_NEG_INFINITY))
					: x(fX)
					, y(fY)
				{
				}
			};

			StreamParticle2D(float fPosX, float fPosY)
				: iIdx(0)
				, iCount(1)
			{
				vPos[0] = Pos(fPosX, fPosY);
			}

			static const uint TRACE_LENGTH = 8;

			Pos vPos[TRACE_LENGTH];
			uint iIdx, iCount;

			inline const Pos& GetCurrentPos() const { return vPos[iIdx]; }
			inline const Pos& GetPosAtIdx(int i) const { return vPos[i]; }
			inline void SetCurrentPos(const Pos& pos)
			{
				iIdx++;
				iIdx %= TRACE_LENGTH;
				if(iCount < TRACE_LENGTH)
					iCount++;

				vPos[iIdx] = pos;
			}
			inline bool OutOfBounds(int iBoundX, int iBoundY) const
			{
				int iTailIdx = iIdx - iCount + 1;
				if(iTailIdx < 0)
					iTailIdx += TRACE_LENGTH;

				float fPosX = vPos[iTailIdx].x;
				float fPosY = vPos[iTailIdx].y;

				return (fPosX < 2 || fPosX > iBoundX - 3 || fPosY < 2 || fPosY > iBoundY - 3);
			}
		};

		class VectorFieldVisualizer2D
		{
		private:
			int miDimX, miDimY;
			RandomGen mRandom;

			vector<StreamParticle2D> maParticles;

		public:
			VectorFieldVisualizer2D()
				: miDimX(0), miDimY(0)
			{
			}
			~VectorFieldVisualizer2D()
			{
			}

			void Initialize(int iDimX, int iDimY, int iParticleCount);
			void Advect(const float* pfVecX, const float* pfVecY, int iDimX, int iDimY, float fDeltaT);
			void AdvectRK4(const float* pfVecX, const float* pfVecY, int iDimX, int iDimY, float fDeltaT);
			void Render(float fScaling, const float* pfColor) const;
			void Release()
			{
				maParticles.clear();
			}
		};

		// 3D Visualizer
		struct StreamParticle3D
		{
			struct Pos3D
			{
				float x;
				float y;
				float z;

				Pos3D(float fX = float(Math::EDX_NEG_INFINITY), float fY = float(Math::EDX_NEG_INFINITY), float fZ = float(Math::EDX_NEG_INFINITY))
					: x(fX)
					, y(fY)
					, z(fZ)
				{
				}
			};

			StreamParticle3D(float fPosX, float fPosY, float fPosZ)
				: iIdx(0)
				, iCount(1)
			{
				vPos[0] = Pos3D(fPosX, fPosY, fPosZ);
			}

			static const uint TRACE_LENGTH = 35;

			Pos3D vPos[TRACE_LENGTH];
			uint iIdx, iCount;

			inline const Pos3D& GetCurrentPos() const { return vPos[iIdx]; }
			inline const Pos3D& GetPosAtIdx(int i) const { return vPos[i]; }
			inline void SetCurrentPos(const Pos3D& pos)
			{
				iIdx++;
				iIdx %= TRACE_LENGTH;
				if(iCount < TRACE_LENGTH)
					iCount++;

				vPos[iIdx] = pos;
			}

			inline bool OutOfBounds(int iBoundX, int iBoundY, int iBoundZ) const
			{
				int iTailIdx = iIdx - iCount + 1;
				if(iTailIdx < 0)
					iTailIdx += TRACE_LENGTH;

				float fPosX = vPos[iTailIdx].x;
				float fPosY = vPos[iTailIdx].y;
				float fPosZ = vPos[iTailIdx].z;

				return (fPosX < 2 || fPosX > iBoundX - 3 || fPosY < 2 || fPosY > iBoundY - 3 || fPosZ < 2 || fPosZ > iBoundZ - 3);
			}
		};

		class VectorFieldVisualizer3D
		{
		private:
			int miDimX, miDimY, miDimZ;

			bool mbPauseAdvect;
			BoundingBox mBBox;
			vector<StreamParticle3D> maParticles;

			RandomGen mRandom;

		public:
			VectorFieldVisualizer3D()
				: miDimX(0), miDimY(0), miDimZ(0)
				, mbPauseAdvect(false)
			{
			}
			~VectorFieldVisualizer3D()
			{
			}

			void Initialize(int iDimX, int iDimY, int iDimZ, int iParticleCount);
			void Advect(const float* pfVecX, const float* pfVecY, const float* pfVecZ, int iDimX, int iDimY, int iDimZ, float fDeltaT);
			//void AdvectRK4(const float* pfVecX, const float* pfVecY, int iDimX, int iDimY, float fDeltaT);
			void Render() const;
			void RenderBBox() const;
			void Release()
			{
				maParticles.clear();
			}

			void TogglePause() { mbPauseAdvect = !mbPauseAdvect; }

		private:
			int Offset(int iX, int iY, int iZ) const { return iZ * miDimX * miDimY + iY * miDimX + iX; }
		};
		
		template<uint Dimension>
		class ParticleVisualizer
		{
		public:
			static void Render(const vector<Particle<Dimension>>& particles, const float fScaling);
		};
	}
}

