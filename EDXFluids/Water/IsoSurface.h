#pragma once

#include "EDXPrerequisites.h"
#include "../Base/ForwardDecl.h"
#include "Memory/Array.h"
#include "Math/Vector.h"

namespace EDX
{
	namespace FluidSim
	{
		struct MCNode
		{
			float		Value;
			Vector3		Pos;
			Vector3		Normal;
		};

		struct MCCell
		{
			MCNode Nodes[8];
		};

		struct IsoVertex
		{
			Vector3	Position;
			Vector3	Normal;
		};

		class IsoSurface
		{
		public:

		private:
			Vector3 mSurfPos;
			int mWidth, mHeight, mDepth;
			int mGridSize;
			MCCell mCell;

			IsoVertex mvVertices[12];

			IsoVertex* mpVertexMem;
			int mVertexCount;

		public:
			float* mpSurfVal;
			float mThreshold;
			IsoSurface(void);
			~IsoSurface(void);

			void Initialize(int sizeX, int sizeY, int sizeZ, int GridSize);
			void UpdateSurface(float* pVal);
			void SetNodeData(int X, int Y, int Z, MCNode& Node);
			void SetCellData(int X, int Y, int Z, int X1, int Y1, int Z1);
			void LerpVertex(IsoVertex& Vertex, MCNode Node1, MCNode Node2);
			void Polygonize();
			void MarchingCube();
			void GenerateObjMesh() const;
			void Render();
			void Release();
		};
	}
}