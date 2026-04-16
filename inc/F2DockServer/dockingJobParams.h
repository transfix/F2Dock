/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.

  F2Dock is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  F2Dock is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include <string>


enum {  LINUX, 
	PRISM2
};

enum {  RECEPTOR_PQR, 
	LIGAND_PQR,
	RECEPTOR_F2D, 
	LIGAND_F2D,
	RECEPTOR_RAWN, 
	LIGAND_RAWN,
	RECEPTOR_QUAD, 
	LIGAND_QUAD,
	DOCKING_OUT,
	RERANKING_OUT
};

enum { 	INIT,
	SUBMITTED, 
	RUNNING, 
	COMPLETED, 
	AVAILABLE
};

enum { 	DOCKING, 
	RERANKING, 
	F2DGEN, 
	QUADGEN
};


typedef struct
{
	int jobId;
	int jobType;

	bool isReceptorPDBAvailable;
	std::string receptorPDBName;

	bool isReceptorPQRAvailable;
	std::string receptorPQRName;

	bool isReceptorF2dAvailable;
	std::string receptorF2dName;

	bool isReceptorRAWNAvailable;
	std::string receptorRAWNName;

	bool isReceptorQuadAvailable;
	std::string receptorQuadName;

	bool isLigandPDBAvailable;
	std::string ligandPDBName;

	bool isLigandPQRAvailable;
	std::string ligandPQRName;

	bool isLigandF2dAvailable;
	std::string ligandF2dName;

	bool isLigandRAWNAvailable;
	std::string ligandRAWNName;

	bool isLigandQuadAvailable;
	std::string ligandQuadName;

	bool isRMSDAvailable;	
	std::string rmsdFileName;

	std::string receptorXYZName;
	std::string ligandXYZName;

	std::string receptorRAWName;
	std::string ligandRAWName;

	std::string receptorIRAWName;
	std::string ligandIRAWName;

	bool isDockingInputAvailable;
	std::string dockingInputFileName;

	bool isRerankingInputAvailable;
	std::string rerankingInputFileName;

	bool isF2dGenInputAvailable;
	std::string f2dGenInputFileName;

	bool isQuadGenInputAvailable;
	std::string quadGenInputFileName;

	bool isDockingOutputAvailable;
	std::string dockingOutputFileName;

	bool isRerankingOutputAvailable;
	std::string rerankingOutputFileName;

	bool performRerank;
	bool storeIntermediateFiles;

	std::string scriptFileName;

	bool isLocal;
	int platform;

	double f2dDim;
	double quadDim;
}DockingJobParams;
