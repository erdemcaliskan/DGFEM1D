#include "InhomogeousField.h"

#define PRNT_FLDOUT_RANDOM 1

OneIHField::OneIHField()
{
	valsAtVertices = true;
	uniformGrid = true;
	isPeriodic = false;
	xm = 0.0;
	xM = 1.0;
	randVariableType = NULL;
	numVertices = 0;
	numSegments = 0;
}

OneIHField::~OneIHField()
{
	if (randVariableType != NULL)
		delete randVariableType;
}

void OneIHField::Read_InstructionsOnly(istream* inConfigPtr, bool* isPeriodicPtr, double* xMPtr, double* xmPtr)
{
	num_Vals_and_x_Provided = false;
	containsRepeatingEndPeriodicVal = false;

	if (isPeriodicPtr != NULL)
		isPeriodic = *isPeriodicPtr;
	if (xMPtr != NULL)
		xM = *xMPtr;
	if (xmPtr != NULL)
		xm = *xmPtr;
	if (inConfigPtr != NULL)
	{
		string buf;
		READ_NSTRING(*inConfigPtr, buf, buf);
		if (buf != "{")
		{
			if (buf == "}")
				return;
			else
			{
				THROW("istream should start with {");
			}
		}
		READ_NSTRING(*inConfigPtr, buf, buf);
		while (buf != "}")
		{
			if (buf == "valsAtVertices")
			{
				READ_NBOOL(*inConfigPtr, buf, valsAtVertices);
			}
			else if (buf == "uniformGrid")
			{
				READ_NBOOL(*inConfigPtr, buf, uniformGrid);
			}
			else if (buf == "isPeriodic")
			{
				bool isPeriodicTmp;
				READ_NBOOL(*inConfigPtr, buf, isPeriodicTmp);
				if (isPeriodicPtr == NULL)
					isPeriodic = isPeriodicTmp;
			}
			else if (buf == "num_Vals_and_x_Provided")
			{
				READ_NBOOL(*inConfigPtr, buf, num_Vals_and_x_Provided);
			}
			else if (buf == "containsRepeatingEndPeriodicVal")
			{
				READ_NBOOL(*inConfigPtr, buf, containsRepeatingEndPeriodicVal);
			}
			else if (buf == "xm")
			{
				double xmTmp;
				READ_NDOUBLE(*inConfigPtr, buf, xmTmp);
				if (xmPtr == NULL)
					xm = xmTmp;
			}
			else if (buf == "xM")
			{
				double xMTmp;
				READ_NDOUBLE(*inConfigPtr, buf, xMTmp);
				if (xMPtr == NULL)
					xM = xMTmp;
			}
			else if (buf == "randVariableType")
			{
				randVariableType = ReadRandomVariableFromFile(*inConfigPtr);
			}
			else
			{
				cout << "buf:\t" << buf << '\n';
				THROW("invalid option\n");
			}
			READ_NSTRING(*inConfigPtr, buf, buf);
		}
	}
}

void OneIHField::Read_DataOnly(istream& inData)
{
	Read_Vals_xs(inData, num_Vals_and_x_Provided, containsRepeatingEndPeriodicVal);
	Finalize_spatialPositions();
	Finalize_Values();
}

void OneIHField::Read_Initialize_OneIHField(istream& inData, istream* inConfigPtr, bool* isPeriodicPtr, double* xMPtr, double* xmPtr, int resolutionFactor, setStatOp_type sso)
{
	Read_InstructionsOnly(inConfigPtr, isPeriodicPtr, xMPtr, xmPtr);
	Read_DataOnly(inData);

#if PRNT_FLDOUT_RANDOM
	static int cntr = 0;
	string cntrStr;
	ostringstream convert;
	convert << cntr;
	cntrStr = convert.str();
//	string name = g_prefileName + "field" + cntrStr + "_initialResolution.txt";
	string name = "field" + cntrStr + "_initialResolution.txt";
	fstream out(name.c_str(), ios::out);
	for (unsigned int i = 0; i < vals.size(); ++i)
		out << vals[i] << '\n';
	out.close();
#endif

	Modify_Resolution(resolutionFactor, sso);

#if PRNT_FLDOUT_RANDOM
//	name = g_prefileName + "field" + cntrStr + "_finalResolution.txt";
	name = "field" + cntrStr + "_finalResolution.txt";
	fstream out2(name.c_str(), ios::out);
	for (unsigned int i = 0; i < vals.size(); ++i)
		out2 << vals[i] << '\n';
	out2.close();
#endif
}

void OneIHField::Modify_Resolution(int resolutionFactor, setStatOp_type sso)
{
	if ((resolutionFactor == 1) || (resolutionFactor == 0) || (resolutionFactor == -1) || (sso == sso_none))
		return;

	vector<double> vals_BK = vals;
	vector<double> xs_BK = xs;
	unsigned int numVertices_BK = numVertices;
	unsigned int numSegments_BK = numSegments;

	if (resolutionFactor > 1)
	{
		unsigned numSegmentsNew = numSegments / resolutionFactor;
		if (numSegmentsNew == 0)
			numSegmentsNew = 1;
		unsigned int numOldSegInNewSeg = numSegments / numSegmentsNew;
		vector<int> subsegmentSizes(numSegmentsNew), starPos(numSegmentsNew + 1);
		vector<vector<double> > weights(numSegmentsNew);
		starPos[0] = 0;
		for (unsigned int i = 0; i < numSegmentsNew - 1; ++i)
		{
			subsegmentSizes[i] = numOldSegInNewSeg;
			starPos[i + 1] = starPos[i] + numOldSegInNewSeg;
		}
		starPos[numSegmentsNew] = numSegments;
		subsegmentSizes[numSegmentsNew - 1] = starPos[numSegmentsNew] - starPos[numSegmentsNew - 1];	//		numSegments - numOldSegInNewSeg * (numSegmentsNew - 1);
		if (!uniformGrid)
			for (unsigned int i = 0; i < numSegmentsNew; ++i)
			{
				unsigned int st = starPos[i], en = starPos[i + 1], szz = en - st;
				vector<double> wts(szz);
				for (unsigned int j = st; j < en; ++j)
					wts[j - st] = xs[j + 1] - xs[j];
				weights.push_back(wts);
			}

		numSegments = numSegmentsNew;
		numVertices = numSegments + 1;
		xs.resize(numVertices);
		for (unsigned int i = 0; i < numVertices; ++i)
			xs[i] = xs_BK[starPos[i]];

		if (valsAtVertices)
			vals.resize(numVertices);
		else
			vals.resize(numSegments);
		for (unsigned int i = 0; i < numSegments; ++i)
		{
			unsigned int st = starPos[i], en = starPos[i + 1], szz = en - st;
			vector<double> valsTmp(szz);
			for (unsigned int j = st; j < en; ++j)
				valsTmp[j - st] = vals_BK[j];
			vals[i] = getWoWOWeightStatvalue(valsTmp, !uniformGrid, weights[i], sso);
		}
	}
	else
	{
		resolutionFactor *= -1;
		numSegments = numSegments_BK * resolutionFactor;
		numVertices = numSegments + 1;
		xs.resize(numVertices);
		if (valsAtVertices)
			vals.resize(numVertices);
		else
			vals.resize(numSegments);
		unsigned int cntr = 0;

		vector<double> facts0(resolutionFactor), facts1(resolutionFactor);
		double fact0, fact1;
		for (unsigned int j = 0; j < (unsigned int)resolutionFactor; ++j)
		{
			facts1[j] = (double)j / (double)resolutionFactor;
			facts0[j] = 1.0 - facts1[j];
		}

		for (unsigned int i = 0; i < numSegments_BK; ++i)
		{
			for (unsigned int j = 0; j < (unsigned int)resolutionFactor; ++j)
			{
				fact0 = facts0[j], fact1 = facts1[j];
				xs[cntr] = xs_BK[i] * fact0 + xs_BK[i + 1] * fact1;
				if (valsAtVertices)
					vals[cntr] = vals_BK[i] * fact0 + vals_BK[i + 1] * fact1;
				else
					vals[cntr] = vals_BK[i];
				++cntr;
			}
		}
		xs[numSegments] = xs_BK[numSegments_BK];
	}
	if (valsAtVertices)
	{
		if (isPeriodic)
			vals[numSegments] = vals[0];
		else
			vals[numSegments] = vals_BK[numSegments_BK];
	}
}

double OneIHField::getVertexValueByIndex(unsigned int index) const
{
	if (valsAtVertices)
		return vals[index];
	unsigned sz = vals.size();
	if (index >= sz)
		index = sz - 1;
	return vals[index];
}

double OneIHField::getSegmentValueByIndex(unsigned int index) const
{
	if (!valsAtVertices)
		return vals[index];
	if (uniformGrid)
		return 0.5 * (vals[index] + vals[index + 1]);
	double x = 0.5 * (xs[index] + xs[index + 1]);
	return getValueByCoord(x);
}

double OneIHField::getValueByCoord(double x) const
{
	const double tol = 1e-12;
	if (uniformGrid)
	{
		double rel;
		unsigned int relInt;
		rel = (x - xm) / (xM - xm) * numSegments;
		bool isBegin = false, isEnd = false;

		if (rel < tol)
		{
			rel = 0.0;
			relInt = 0;
			return vals[0];
		}
		else if (numSegments - rel < tol)
		{
			rel = 1.0;
			relInt = numSegments;
			return vals[vals.size() - 1];
		}
		relInt = (int)floor(rel);
		if (!valsAtVertices)
			return vals[relInt];
		double vm = vals[relInt];
		double vM = vals[relInt + 1];
		double fM = rel - relInt, fm = 1.0 - fM;
		return vm * fm + vM * fM;
	}
	// non-uniform grid
	double relTol = tol * (xM - xm);
	int sz_vl = vals.size();
	if (xM - x < relTol)
		return vals[sz_vl - 1];
	if (x - xm < relTol)
		return vals[0];
	int i;
	for (i = 0; i < sz_vl; ++i)
	{
		if (x < xs[i])
			break;
	}
	if (!valsAtVertices)
		return vals[i - 1];
	double vm = vals[i - 1], vM = vals[i];
	double fM = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
	double fm = 1.0 - fM;
	return fm * vm + fM * vM;
}

void OneIHField::Output(ostream& out)
{
	unsigned int sz_vals = vals.size();
	out << sz_vals << '\n';
	for (unsigned int i = 0; i < sz_vals; ++i)
		out << vals[i] << '\t';
	out << '\n';

	unsigned int sz_xs = xs.size();
	out << sz_xs << '\n';
	for (unsigned int i = 0; i < sz_xs; ++i)
		out << xs[i] << '\t';
	out << '\n';

	out << "numVertices\t" << numVertices << '\n';
	out << "numSegments\t" << numSegments << '\n';
	out << "valsAtVertices\t" << valsAtVertices << '\n';
	out << "uniformGrid\t" << uniformGrid << '\n';
	out << "isPeriodic\t" << isPeriodic << '\n';
	out << "xm\t" << xm << '\n';
	out << "xM\t" << xM << '\n';
	if (randVariableType != NULL)
	{
		randVariableType->paras.WriteData(out);
		out << '\n';
	}
}

OneIHField::OneIHField(const OneIHField& other)
{
	randVariableType = NULL;
	(*this) = other;
}

OneIHField& OneIHField::operator=(const OneIHField& other)
{
	valsAtVertices = other.valsAtVertices;
	uniformGrid = other.uniformGrid;
	isPeriodic = other.isPeriodic;
	xm = other.xm;
	xM = other.xM;
	randVariableType = other.randVariableType->CreateCopy();
	vals = other.vals;
	xs = other.xs;
	numVertices = other.numVertices;
	numSegments = other.numSegments;
	return *this;
}

void OneIHField::Read_Vals_xs(istream & in, bool num_vals_and_x_Provided, bool containsRepeatingEndPeriodicVal)
{
	bool addLastPeriodicValue = isPeriodic && valsAtVertices && !containsRepeatingEndPeriodicVal;
	if (num_vals_and_x_Provided)
	{
		if (uniformGrid)
		{
			int num_vals;
			in >> num_vals;

			vals.resize(num_vals);
			for (int i = 0; i < num_vals; ++i)
				in >> vals[i];
		}
		else
		{
			int num_xs;
			in >> num_xs;

			xs.resize(num_xs);
			for (int i = 0; i < num_xs; ++i)
				in >> xs[i];

			int num_vals;
			in >> num_vals;
			vals.resize(num_vals);
			for (int i = 0; i < num_vals; ++i)
				in >> vals[i];
		}
	}
	else
	{
		double tmp;
		vector<double> tmps;
		in >> tmp;
		while (!in.eof())
		{
			tmps.push_back(tmp);
			in >> tmp;
		}
		unsigned int sz = tmps.size();

		if (uniformGrid)
			vals = tmps;
		else
		{
			unsigned int num_xs;
			if (!valsAtVertices || addLastPeriodicValue) // 2nV - 1 values read (nV number of vertices)
				num_xs = (sz + 1) / 2;
			else // 2nV values read (nV number of vertices)
				num_xs = sz / 2;

			xs.resize(num_xs);
			for (unsigned int i = 0; i < num_xs; ++i)
				xs[i] = tmps[i];
			unsigned int num_vals = sz - num_xs;
			vals.resize(num_vals);
			for (unsigned int i = 0; i < num_vals; ++i)
				vals[i] = tmps[i + num_xs];
		}
	}
	if (addLastPeriodicValue)
		vals.push_back(vals[0]);
	unsigned int num_vals = vals.size();

	if (valsAtVertices)
	{
		numVertices = num_vals;
		numSegments = numVertices - 1;
	}
	else
	{
		numSegments = num_vals;
		numVertices = numSegments + 1;
	}
}

void OneIHField::Finalize_spatialPositions()
{
	if (uniformGrid)
	{
		if (xs.size() != numVertices)
			xs.resize(numVertices);

		double delV = (xM - xm) / numSegments;
		for (unsigned int i = 0; i < numVertices; ++i)
			xs[i] = xm + delV * i;
	}
	else
	{
		if (xs.size() != numVertices)
		{
			cout << "xs size\t" << xs.size();
			THROW("Invalid size\t");
		}
		double xmxs = xs[0], xMxs = xs[numVertices - 1];
		double delx = (xM - xm);
		double b = delx / (xMxs - xmxs);
		double a = xm - b * xmxs;
		for (unsigned int i = 0; i < numVertices; ++i)
			xs[i] = a + b * xs[i];
	}
}

void OneIHField::Finalize_Values()
{
	if (randVariableType == NULL)
		return;
	for (unsigned int i = 0; i < vals.size(); ++i)
		vals[i] = randVariableType->TurnStandardNormalValue2ThisRandom(vals[i]);
}

void TestInhomogeneousField(string baseNameWOExt)
{
	string fileNameDat = baseNameWOExt + ".txt";
	string fileNameInstructions = baseNameWOExt + ".inst";
	fstream inData(fileNameDat.c_str(), ios::in);
	fstream inConfig(fileNameInstructions.c_str(), ios::in);
	OneIHField oihf;
	bool* isPeriodicPtr = NULL;
	double* xMPtr = NULL; double* xmPtr = NULL;
	oihf.Read_Initialize_OneIHField(inData, &inConfig, isPeriodicPtr, xMPtr, xmPtr, 1, sso_mean_arithmetic);

	string fileNameOut = baseNameWOExt + ".out";
	fstream out(fileNameOut.c_str(), ios::out);
	oihf.Output(out);

	string fileNameVals = baseNameWOExt + "_vals.txt";
	fstream outv(fileNameVals.c_str(), ios::out);
	unsigned int numVals = oihf.getNumValues();
	for (unsigned int i = 0; i < numVals; ++i)
		outv << oihf.getVertexValueByIndex(i) << '\t';
	outv << '\n';
	out << "numVals\t" << numVals << '\n';

	unsigned int numSegments = oihf.getNumSegments();
	double xm, xM;
	oihf.get_domain_range(xm, xM);
	unsigned int facti = 100;
	double del = (xM - xm) / numSegments / facti;
	unsigned numPts = facti * numSegments + 1;
//	vector<double> xss(numPts), vls(numPts);
	double x, v;

	string fileNameValsMatlab = baseNameWOExt + "_x_v.txt";
	fstream outxv(fileNameValsMatlab.c_str(), ios::out);

	for (unsigned int i = 0; i < numPts; ++i)
	{
		x = xm + i * del;
		v = oihf.getValueByCoord(x);
		outxv << x << '\t' << v << '\n';
	}
}
