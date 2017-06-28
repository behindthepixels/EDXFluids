#include "VectorVisualizer.h"
#include "Math/EDXMath.h"
#include "Math/BoundingBox.h"
#include "Core/Random.h"

#include "Graphics/OpenGL.h"

#include <ppl.h>
using namespace Concurrency;

namespace EDX
{
	namespace FluidSim
	{
		void VectorFieldVisualizer2D::Initialize(int iDimX, int iDimY, int iParticleCount)
		{
			miDimX = iDimX;
			miDimY = iDimY;

			RandomGen rng;
			for(int i = 0; i < iParticleCount; i++)
			{
				maParticles.push_back(StreamParticle2D(rng.Float() * miDimX, rng.Float() * miDimY));
			}
		}

		void VectorFieldVisualizer2D::Advect(const float* pfVecX, const float* pfVecY, int iDimX, int iDimY, float fDeltaT)
		{
			assert(miDimX == iDimX && miDimY == iDimY);

			float fDT = fDeltaT * miDimX;
			for(int i = 0; i < maParticles.size(); i++)
			{
				StreamParticle2D& particle = maParticles[i];

				StreamParticle2D::Pos pos = particle.GetCurrentPos();

				float fPosX = Math::Clamp(pos.x, 0.0f, float(miDimX - 1));
				int iX0 = Math::FloorToInt(fPosX);
				int iX1 = Math::Min(iX0 + 1, miDimX - 1);
				int iXm1 = Math::Max(iX0 - 1, 0);
				int iX2 = Math::Min(iX1 + 1, miDimX - 1);

				float fPosY = Math::Clamp(pos.y, 0.0f, float(miDimY - 1));
				int iY0 = Math::FloorToInt(fPosY);
				int iY1 = Math::Min(iY0 + 1, miDimY - 1);
				int iYm1 = Math::Max(iY0 - 1, 0);
				int iY2 = Math::Min(iY1 + 1, miDimY - 1);

				float fLinX = Math::LinStep(fPosX, iX0, iX1);
				float fLinY = Math::LinStep(fPosY, iY0, iY1);

				float fTempX0 = Math::Lerp(pfVecX[iY0 * miDimX + iX0], pfVecX[iY0 * miDimX + iX1], fLinX);
				float fTempX1 = Math::Lerp(pfVecX[iY1 * miDimX + iX0], pfVecX[iY1 * miDimX + iX1], fLinX);
				float fVelX = Math::Lerp(fTempX0, fTempX1, fLinY);

				float fTempY0 = Math::Lerp(pfVecY[iY0 * miDimX + iX0], pfVecY[iY0 * miDimX + iX1], fLinX);
				float fTempY1 = Math::Lerp(pfVecY[iY1 * miDimX + iX0], pfVecY[iY1 * miDimX + iX1], fLinX);
				float fVelY = Math::Lerp(fTempY0, fTempY1, fLinY);

				fVelX = Math::Clamp(fVelX, -1.0f, 1.0f);
				fVelY = Math::Clamp(fVelY, -1.0f, 1.0f);

				float fNewX = fPosX + fVelX * fDT;
				float fNewY = fPosY + fVelY * fDT;

				fNewX = Math::Clamp(fNewX, 0, miDimX - 1);
				fNewY = Math::Clamp(fNewY, 0, miDimY - 1);

				particle.SetCurrentPos(StreamParticle2D::Pos(fNewX, fNewY));
			}
		}

		void VectorFieldVisualizer2D::AdvectRK4(const float* pfVecX, const float* pfVecY, int iDimX, int iDimY, float fDeltaT)
		{
			//assert(miDimX == iDimX && miDimY == iDimY);

			float fDT = fDeltaT * iDimY;
			float fDT_2 = fDT * 0.5f;
			for(int i = 0; i < maParticles.size(); i++)
			{
				StreamParticle2D& particle = maParticles[i];

				StreamParticle2D::Pos pos = particle.GetCurrentPos();

				float fPosX = Math::Clamp(pos.x, 0.0f, float(miDimX - 1));
				int iX0 = Math::FloorToInt(fPosX);
				int iX1 = Math::Min(iX0 + 1, miDimX - 1);

				float fPosY = Math::Clamp(pos.y, 0.0f, float(miDimY - 1));
				int iY0 = Math::FloorToInt(fPosY);
				int iY1 = Math::Min(iY0 + 1, miDimY - 1);

				float fLinX = Math::LinStep(fPosX, iX0, iX1);
				float fLinY = Math::LinStep(fPosY, iY0, iY1);

				float fTempX0 = Math::Lerp(pfVecX[iY0 * miDimX + iX0], pfVecX[iY0 * miDimX + iX1], fLinX);
				float fTempX1 = Math::Lerp(pfVecX[iY1 * miDimX + iX0], pfVecX[iY1 * miDimX + iX1], fLinX);
				float fVelX = Math::Lerp(fTempX0, fTempX1, fLinY);

				float fTempY0 = Math::Lerp(pfVecY[iY0 * miDimX + iX0], pfVecY[iY0 * miDimX + iX1], fLinX);
				float fTempY1 = Math::Lerp(pfVecY[iY1 * miDimX + iX0], pfVecY[iY1 * miDimX + iX1], fLinX);
				float fVelY = Math::Lerp(fTempY0, fTempY1, fLinY);

				// 1st
				float fPosX1 = fPosX + fVelX * fDT_2;
				float fPosY1 = fPosY + fVelY * fDT_2;

				fPosX1 = Math::Clamp(fPosX1, 0, miDimX - 1);
				fPosY1 = Math::Clamp(fPosY1, 0, miDimY - 1);

				iX0 = Math::FloorToInt(fPosX1);
				iX1 = Math::Min(iX0 + 1, miDimX - 1);

				iY0 = Math::FloorToInt(fPosY1);
				iY1 = Math::Min(iY0 + 1, miDimY - 1);

				fLinX = Math::LinStep(fPosX1, iX0, iX1);
				fLinY = Math::LinStep(fPosY1, iY0, iY1);

				fTempX0 = Math::Lerp(pfVecX[iY0 * miDimX + iX0], pfVecX[iY0 * miDimX + iX1], fLinX);
				fTempX1 = Math::Lerp(pfVecX[iY1 * miDimX + iX0], pfVecX[iY1 * miDimX + iX1], fLinX);
				float fVelX1 = Math::Lerp(fTempX0, fTempX1, fLinY);

				fTempY0 = Math::Lerp(pfVecY[iY0 * miDimX + iX0], pfVecY[iY0 * miDimX + iX1], fLinX);
				fTempY1 = Math::Lerp(pfVecY[iY1 * miDimX + iX0], pfVecY[iY1 * miDimX + iX1], fLinX);
				float fVelY1 = Math::Lerp(fTempY0, fTempY1, fLinY);

				// 2nd
				float fPosX2 = fPosX + fVelX1 * fDT_2;
				float fPosY2 = fPosY + fVelY1 * fDT_2;

				fPosX2 = Math::Clamp(fPosX2, 0, miDimX - 1);
				fPosY2 = Math::Clamp(fPosY2, 0, miDimY - 1);

				iX0 = Math::FloorToInt(fPosX2);
				iX1 = Math::Min(iX0 + 1, miDimX - 1);

				iY0 = Math::FloorToInt(fPosY2);
				iY1 = Math::Min(iY0 + 1, miDimY - 1);

				fLinX = Math::LinStep(fPosX2, iX0, iX1);
				fLinY = Math::LinStep(fPosY2, iY0, iY1);

				fTempX0 = Math::Lerp(pfVecX[iY0 * miDimX + iX0], pfVecX[iY0 * miDimX + iX1], fLinX);
				fTempX1 = Math::Lerp(pfVecX[iY1 * miDimX + iX0], pfVecX[iY1 * miDimX + iX1], fLinX);
				float fVelX2 = Math::Lerp(fTempX0, fTempX1, fLinY);

				fTempY0 = Math::Lerp(pfVecY[iY0 * miDimX + iX0], pfVecY[iY0 * miDimX + iX1], fLinX);
				fTempY1 = Math::Lerp(pfVecY[iY1 * miDimX + iX0], pfVecY[iY1 * miDimX + iX1], fLinX);
				float fVelY2 = Math::Lerp(fTempY0, fTempY1, fLinY);

				// 3rd
				float fPosX3 = fPosX + fVelX2 * fDT;
				float fPosY3 = fPosY + fVelY2 * fDT;

				fPosX3 = Math::Clamp(fPosX3, 0, miDimX - 1);
				fPosY3 = Math::Clamp(fPosY3, 0, miDimY - 1);

				iX0 = Math::FloorToInt(fPosX3);
				iX1 = Math::Min(iX0 + 1, miDimX - 1);

				iY0 = Math::FloorToInt(fPosY3);
				iY1 = Math::Min(iY0 + 1, miDimY - 1);

				fLinX = Math::LinStep(fPosX3, iX0, iX1);
				fLinY = Math::LinStep(fPosY3, iY0, iY1);

				fTempX0 = Math::Lerp(pfVecX[iY0 * miDimX + iX0], pfVecX[iY0 * miDimX + iX1], fLinX);
				fTempX1 = Math::Lerp(pfVecX[iY1 * miDimX + iX0], pfVecX[iY1 * miDimX + iX1], fLinX);
				float fVelX3 = Math::Lerp(fTempX0, fTempX1, fLinY);

				fTempY0 = Math::Lerp(pfVecY[iY0 * miDimX + iX0], pfVecY[iY0 * miDimX + iX1], fLinX);
				fTempY1 = Math::Lerp(pfVecY[iY1 * miDimX + iX0], pfVecY[iY1 * miDimX + iX1], fLinX);
				float fVelY3 = Math::Lerp(fTempY0, fTempY1, fLinY);

				// 4th
				float fPosX4 = fPosX + (fVelX + 2.0f * fVelX1 + 2.0f * fVelX2 + fVelX3) / 6.0f * fDT;
				float fPosY4 = fPosY + (fVelY + 2.0f * fVelY1 + 2.0f * fVelY2 + fVelY3) / 6.0f * fDT;

				fPosX4 = Math::Clamp(fPosX4, 0, miDimX - 1);
				fPosY4 = Math::Clamp(fPosY4, 0, miDimY - 1);

				if(particle.OutOfBounds(miDimX, miDimY))
				{
					fPosX4 = mRandom.Float() * miDimX;
					fPosY4 = mRandom.Float() * miDimY;

					particle.iCount = 0;
				}

				particle.SetCurrentPos(StreamParticle2D::Pos(fPosX4, fPosY4));
			}
		}

		void VectorFieldVisualizer2D::Render(float fScaling, const float* pfColor) const
		{
			for(int i = 0; i < maParticles.size(); i++)
			{
				const StreamParticle2D& particle = maParticles[i];

				int iIndex = particle.iIdx;

				glColor3f(0.0f, 0.0f, 0.0f);
				glBegin(GL_LINE_STRIP);
				for(auto j = 0; j < particle.iCount - 1; j++)
				{
					const StreamParticle2D::Pos& pos = particle.GetPosAtIdx(iIndex);

					int iX = Math::RoundToInt(pos.x);
					int iY = Math::RoundToInt(pos.y);
					if(!pfColor || pfColor[iY * miDimX + iX] <= 0.0f)
					{
						glColor3f(0.0f, 0.0f, 1.0f);
					}
					else
					{
						glColor3f(1.0f, 0.0f, 0.0f);
					}

					glVertex2f((pos.x + 0.5f) * fScaling, (pos.y + 0.5f) * fScaling);

					iIndex--;
					if(iIndex < 0)
					{
						iIndex = StreamParticle2D::TRACE_LENGTH - 1;
					}
				}
				glEnd();
			}
		}

		// -------------------------------------------------------------------------------------
		// 3D particle trace visualizer
		// -------------------------------------------------------------------------------------
		void VectorFieldVisualizer3D::Initialize(int iDimX, int iDimY, int iDimZ, int iParticleCount)
		{
			miDimX = iDimX;
			miDimY = iDimY;
			miDimZ = iDimZ;

			//mBBox = BoundingBox(Vector3::ZERO, Vector3(miDimX, miDimY, miDimZ));
			mBBox.mMin = Vector3::ZERO;
			mBBox.mMax = Vector3(miDimX, miDimY, miDimZ);

			RandomGen rng;
			for(int i = 0; i < iParticleCount; i++)
			{
				maParticles.push_back(StreamParticle3D(rng.Float() * miDimX, rng.Float() * miDimY, rng.Float() * miDimZ));
			}
		}
		void VectorFieldVisualizer3D::Advect(const float* pfVecX, const float* pfVecY, const float* pfVecZ, int iDimX, int iDimY, int iDimZ, float fDeltaT)
		{
			if(mbPauseAdvect)				
			{
				return;
			}

			float fDT = fDeltaT * iDimY;
			float fDT_2 = fDT * 0.5f;
			int iSize = maParticles.size();
			parallel_for(0, iSize - 1, [&](int i)
			//for(int i = 0; i < maParticles.size(); i++)
			{
				StreamParticle3D& particle = maParticles[i];

				StreamParticle3D::Pos3D pos = particle.GetCurrentPos();

				float fPosX = Math::Clamp(pos.x, 1.0f, float(miDimX - 1));
				int iX0 = Math::FloorToInt(fPosX);
				int iX1 = Math::Min(iX0 + 1, miDimX - 2);

				float fPosY = Math::Clamp(pos.y, 1.0f, float(miDimY - 1));
				int iY0 = Math::FloorToInt(fPosY);
				int iY1 = Math::Min(iY0 + 1, miDimY - 2);

				float fPosZ = Math::Clamp(pos.z, 1.0f, float(miDimZ - 1));
				int iZ0 = Math::FloorToInt(fPosZ);
				int iZ1 = Math::Min(iZ0 + 1, miDimZ - 2);

				float fLinX = Math::LinStep(fPosX, iX0, iX1);
				float fLinY = Math::LinStep(fPosY, iY0, iY1);
				float fLinZ = Math::LinStep(fPosZ, iZ0, iZ1);

				// X
				float fTempX00 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ0)], pfVecX[Offset(iX1, iY0, iZ0)], fLinX);
				float fTempX01 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ1)], pfVecX[Offset(iX1, iY0, iZ1)], fLinX);
				float fTempX10 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ0)], pfVecX[Offset(iX1, iY1, iZ0)], fLinX);
				float fTempX11 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ1)], pfVecX[Offset(iX1, iY1, iZ1)], fLinX);

				float fTempX0 = Math::Lerp(fTempX00, fTempX10, fLinY);
				float fTempX1 = Math::Lerp(fTempX01, fTempX11, fLinY);

				float fVelX = Math::Lerp(fTempX0, fTempX1, fLinZ);

				// Y
				float fTempY00 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ0)], pfVecY[Offset(iX1, iY0, iZ0)], fLinX);
				float fTempY01 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ1)], pfVecY[Offset(iX1, iY0, iZ1)], fLinX);
				float fTempY10 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ0)], pfVecY[Offset(iX1, iY1, iZ0)], fLinX);
				float fTempY11 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ1)], pfVecY[Offset(iX1, iY1, iZ1)], fLinX);

				float fTempY0 = Math::Lerp(fTempY00, fTempY10, fLinY);
				float fTempY1 = Math::Lerp(fTempY01, fTempY11, fLinY);

				float fVelY = Math::Lerp(fTempY0, fTempY1, fLinZ);

				// Z
				float fTempZ00 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ0)], pfVecZ[Offset(iX1, iY0, iZ0)], fLinX);
				float fTempZ01 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ1)], pfVecZ[Offset(iX1, iY0, iZ1)], fLinX);
				float fTempZ10 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ0)], pfVecZ[Offset(iX1, iY1, iZ0)], fLinX);
				float fTempZ11 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ1)], pfVecZ[Offset(iX1, iY1, iZ1)], fLinX);

				float fTempZ0 = Math::Lerp(fTempZ00, fTempZ10, fLinY);
				float fTempZ1 = Math::Lerp(fTempZ01, fTempZ11, fLinY);

				float fVelZ = Math::Lerp(fTempZ0, fTempZ1, fLinZ);

				// 1st
				float fPosX1 = fPosX + fVelX * fDT_2;
				float fPosY1 = fPosY + fVelY * fDT_2;
				float fPosZ1 = fPosZ + fVelZ * fDT_2;

				fPosX1 = Math::Clamp(fPosX1, 1.0f, float(miDimX - 1));
				fPosY1 = Math::Clamp(fPosY1, 1.0f, float(miDimY - 1));
				fPosZ1 = Math::Clamp(fPosZ1, 1.0f, float(miDimZ - 1));

				iX0 = Math::FloorToInt(fPosX1);
				iX1 = Math::Min(iX0 + 1, miDimX - 2);

				iY0 = Math::FloorToInt(fPosY1);
				iY1 = Math::Min(iY0 + 1, miDimY - 2);

				iZ0 = Math::FloorToInt(fPosZ1);
				iZ1 = Math::Min(iZ0 + 1, miDimZ - 2);

				fLinX = Math::LinStep(fPosX1, iX0, iX1);
				fLinY = Math::LinStep(fPosY1, iY0, iY1);
				fLinZ = Math::LinStep(fPosZ1, iZ0, iZ1);

				// X
				fTempX00 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ0)], pfVecX[Offset(iX1, iY0, iZ0)], fLinX);
				fTempX01 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ1)], pfVecX[Offset(iX1, iY0, iZ1)], fLinX);
				fTempX10 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ0)], pfVecX[Offset(iX1, iY1, iZ0)], fLinX);
				fTempX11 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ1)], pfVecX[Offset(iX1, iY1, iZ1)], fLinX);

				fTempX0 = Math::Lerp(fTempX00, fTempX10, fLinY);
				fTempX1 = Math::Lerp(fTempX01, fTempX11, fLinY);

				float fVelX1 = Math::Lerp(fTempX0, fTempX1, fLinZ);

				// Y
				fTempY00 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ0)], pfVecY[Offset(iX1, iY0, iZ0)], fLinX);
				fTempY01 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ1)], pfVecY[Offset(iX1, iY0, iZ1)], fLinX);
				fTempY10 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ0)], pfVecY[Offset(iX1, iY1, iZ0)], fLinX);
				fTempY11 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ1)], pfVecY[Offset(iX1, iY1, iZ1)], fLinX);

				fTempY0 = Math::Lerp(fTempY00, fTempY10, fLinY);
				fTempY1 = Math::Lerp(fTempY01, fTempY11, fLinY);

				float fVelY1 = Math::Lerp(fTempY0, fTempY1, fLinZ);

				// Z
				fTempZ00 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ0)], pfVecZ[Offset(iX1, iY0, iZ0)], fLinX);
				fTempZ01 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ1)], pfVecZ[Offset(iX1, iY0, iZ1)], fLinX);
				fTempZ10 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ0)], pfVecZ[Offset(iX1, iY1, iZ0)], fLinX);
				fTempZ11 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ1)], pfVecZ[Offset(iX1, iY1, iZ1)], fLinX);

				fTempZ0 = Math::Lerp(fTempZ00, fTempZ10, fLinY);
				fTempZ1 = Math::Lerp(fTempZ01, fTempZ11, fLinY);

				float fVelZ1 = Math::Lerp(fTempZ0, fTempZ1, fLinZ);


				// 2nd
				float fPosX2 = fPosX + fVelX1 * fDT_2;
				float fPosY2 = fPosY + fVelY1 * fDT_2;
				float fPosZ2 = fPosZ + fVelZ1 * fDT_2;

				fPosX2 = Math::Clamp(fPosX2, 1.0f, float(miDimX - 1));
				fPosY2 = Math::Clamp(fPosY2, 1.0f, float(miDimY - 1));
				fPosZ2 = Math::Clamp(fPosZ2, 1.0f, float(miDimZ - 1));

				iX0 = Math::FloorToInt(fPosX2);
				iX1 = Math::Min(iX0 + 1, miDimX - 2);

				iY0 = Math::FloorToInt(fPosY2);
				iY1 = Math::Min(iY0 + 1, miDimY - 2);

				iZ0 = Math::FloorToInt(fPosZ2);
				iZ1 = Math::Min(iZ0 + 1, miDimZ - 2);

				fLinX = Math::LinStep(fPosX2, iX0, iX1);
				fLinY = Math::LinStep(fPosY2, iY0, iY1);
				fLinZ = Math::LinStep(fPosZ2, iZ0, iZ1);

				// X
				fTempX00 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ0)], pfVecX[Offset(iX1, iY0, iZ0)], fLinX);
				fTempX01 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ1)], pfVecX[Offset(iX1, iY0, iZ1)], fLinX);
				fTempX10 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ0)], pfVecX[Offset(iX1, iY1, iZ0)], fLinX);
				fTempX11 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ1)], pfVecX[Offset(iX1, iY1, iZ1)], fLinX);

				fTempX0 = Math::Lerp(fTempX00, fTempX10, fLinY);
				fTempX1 = Math::Lerp(fTempX01, fTempX11, fLinY);

				float fVelX2 = Math::Lerp(fTempX0, fTempX1, fLinZ);

				// Y
				fTempY00 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ0)], pfVecY[Offset(iX1, iY0, iZ0)], fLinX);
				fTempY01 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ1)], pfVecY[Offset(iX1, iY0, iZ1)], fLinX);
				fTempY10 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ0)], pfVecY[Offset(iX1, iY1, iZ0)], fLinX);
				fTempY11 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ1)], pfVecY[Offset(iX1, iY1, iZ1)], fLinX);

				fTempY0 = Math::Lerp(fTempY00, fTempY10, fLinY);
				fTempY1 = Math::Lerp(fTempY01, fTempY11, fLinY);

				float fVelY2 = Math::Lerp(fTempY0, fTempY1, fLinZ);

				// Z
				fTempZ00 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ0)], pfVecZ[Offset(iX1, iY0, iZ0)], fLinX);
				fTempZ01 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ1)], pfVecZ[Offset(iX1, iY0, iZ1)], fLinX);
				fTempZ10 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ0)], pfVecZ[Offset(iX1, iY1, iZ0)], fLinX);
				fTempZ11 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ1)], pfVecZ[Offset(iX1, iY1, iZ1)], fLinX);

				fTempZ0 = Math::Lerp(fTempZ00, fTempZ10, fLinY);
				fTempZ1 = Math::Lerp(fTempZ01, fTempZ11, fLinY);

				float fVelZ2 = Math::Lerp(fTempZ0, fTempZ1, fLinZ);

				// 3rd
				float fPosX3 = fPosX + fVelX2 * fDT;
				float fPosY3 = fPosY + fVelY2 * fDT;
				float fPosZ3 = fPosZ + fVelZ2 * fDT;

				fPosX3 = Math::Clamp(fPosX3, 1.0f, float(miDimX - 1));
				fPosY3 = Math::Clamp(fPosY3, 1.0f, float(miDimY - 1));
				fPosZ3 = Math::Clamp(fPosZ3, 1.0f, float(miDimZ - 1));

				iX0 = Math::FloorToInt(fPosX3);
				iX1 = Math::Min(iX0 + 1, miDimX - 2);

				iY0 = Math::FloorToInt(fPosY3);
				iY1 = Math::Min(iY0 + 1, miDimY - 2);

				iZ0 = Math::FloorToInt(fPosZ3);
				iZ1 = Math::Min(iZ0 + 1, miDimZ - 2);

				fLinX = Math::LinStep(fPosX3, iX0, iX1);
				fLinY = Math::LinStep(fPosY3, iY0, iY1);
				fLinZ = Math::LinStep(fPosZ3, iZ0, iZ1);

				// X
				fTempX00 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ0)], pfVecX[Offset(iX1, iY0, iZ0)], fLinX);
				fTempX01 = Math::Lerp(pfVecX[Offset(iX0, iY0, iZ1)], pfVecX[Offset(iX1, iY0, iZ1)], fLinX);
				fTempX10 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ0)], pfVecX[Offset(iX1, iY1, iZ0)], fLinX);
				fTempX11 = Math::Lerp(pfVecX[Offset(iX0, iY1, iZ1)], pfVecX[Offset(iX1, iY1, iZ1)], fLinX);

				fTempX0 = Math::Lerp(fTempX00, fTempX10, fLinY);
				fTempX1 = Math::Lerp(fTempX01, fTempX11, fLinY);

				float fVelX3 = Math::Lerp(fTempX0, fTempX1, fLinZ);

				// Y
				fTempY00 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ0)], pfVecY[Offset(iX1, iY0, iZ0)], fLinX);
				fTempY01 = Math::Lerp(pfVecY[Offset(iX0, iY0, iZ1)], pfVecY[Offset(iX1, iY0, iZ1)], fLinX);
				fTempY10 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ0)], pfVecY[Offset(iX1, iY1, iZ0)], fLinX);
				fTempY11 = Math::Lerp(pfVecY[Offset(iX0, iY1, iZ1)], pfVecY[Offset(iX1, iY1, iZ1)], fLinX);

				fTempY0 = Math::Lerp(fTempY00, fTempY10, fLinY);
				fTempY1 = Math::Lerp(fTempY01, fTempY11, fLinY);

				float fVelY3 = Math::Lerp(fTempY0, fTempY1, fLinZ);

				// Z
				fTempZ00 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ0)], pfVecZ[Offset(iX1, iY0, iZ0)], fLinX);
				fTempZ01 = Math::Lerp(pfVecZ[Offset(iX0, iY0, iZ1)], pfVecZ[Offset(iX1, iY0, iZ1)], fLinX);
				fTempZ10 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ0)], pfVecZ[Offset(iX1, iY1, iZ0)], fLinX);
				fTempZ11 = Math::Lerp(pfVecZ[Offset(iX0, iY1, iZ1)], pfVecZ[Offset(iX1, iY1, iZ1)], fLinX);

				fTempZ0 = Math::Lerp(fTempZ00, fTempZ10, fLinY);
				fTempZ1 = Math::Lerp(fTempZ01, fTempZ11, fLinY);

				float fVelZ3 = Math::Lerp(fTempZ0, fTempZ1, fLinZ);

				// 4th
				float fPosX4 = fPosX + (fVelX + 2.0f * fVelX1 + 2.0f * fVelX2 + fVelX3) / 6.0f * fDT;
				float fPosY4 = fPosY + (fVelY + 2.0f * fVelY1 + 2.0f * fVelY2 + fVelY3) / 6.0f * fDT;
				float fPosZ4 = fPosZ + (fVelZ + 2.0f * fVelZ1 + 2.0f * fVelZ2 + fVelZ3) / 6.0f * fDT;

				if(particle.OutOfBounds(miDimX, miDimY, miDimZ))
				{
					fPosX4 = mRandom.Float() * miDimX;
					fPosY4 = mRandom.Float() * miDimY;
					fPosZ4 = mRandom.Float() * miDimZ;
					
					particle.iCount = 0;
				}

				fPosX4 = Math::Clamp(fPosX4, 0, miDimX - 1);
				fPosY4 = Math::Clamp(fPosY4, 0, miDimY - 1);
				fPosZ4 = Math::Clamp(fPosZ4, 0, miDimZ - 1);

				particle.SetCurrentPos(StreamParticle3D::Pos3D(fPosX4, fPosY4, fPosZ4));
			});
		}

		void VectorFieldVisualizer3D::Render() const
		{
			for(int i = 0; i < maParticles.size(); i++)
			{
				const StreamParticle3D& particle = maParticles[i];

				int iIndex = particle.iIdx;
				if(particle.iCount <= 1)
				{
					continue;
				}

				glColor3f(0.0f, 0.0f, 0.0f);
				glBegin(GL_LINE_STRIP);
				for(int j = 0; j < particle.iCount - 1; j++)
				{
					const StreamParticle3D::Pos3D& pos = particle.GetPosAtIdx(iIndex);
					int iIndexPrev = iIndex - 1;
					if(iIndexPrev < 0)
					{
						iIndexPrev = StreamParticle3D::TRACE_LENGTH - 1;
					}

					const StreamParticle3D::Pos3D& pos2 = particle.GetPosAtIdx(iIndexPrev);
					
					glTexCoord3f(pos.x - pos2.x, pos.y - pos2.y, pos.z - pos2.z);
					glVertex3f(pos.x, pos.y, pos.z);

					iIndex--;
					if(iIndex < 0)
					{
						iIndex = StreamParticle3D::TRACE_LENGTH - 1;
					}
				}
				glEnd();
			}

		}

		void VectorFieldVisualizer3D::RenderBBox() const
		{
			glColor3f(1.0f, 1.0f, 1.0f);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			//mBBox.RenderInGL();
		}
		
		void ParticleVisualizer<2>::Render(const vector<Particle<2>>& particles, const float fScaling)
		{
			glBegin(GL_POINTS);

			glColor3f(0.0f, 0.0f, 0.0f);
			for(const auto& particle : particles)
				glVertex2f((particle.vPos.x + 0.5f) * fScaling, (particle.vPos.y + 0.5f) * fScaling);

			glEnd();
		}
		
		void ParticleVisualizer<3>::Render(const vector<Particle<3>>& particles, const float fScaling)
		{
			glBegin(GL_POINTS);
			
			glColor3f(1.0f, 1.0f, 1.0f);
			for(const auto& particle : particles)
				glVertex3f((particle.vPos.x + 0.5f) * fScaling, (particle.vPos.y + 0.5f) * fScaling, (particle.vPos.z + 0.5f) * fScaling);

			glEnd();
		}
	}
}