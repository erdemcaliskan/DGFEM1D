#include "RandomVariable.h"

string getName(prob_distrib_type dat)
{
	if (dat == pUniform)
		return "Uniform";
	if (dat == pNormal)
		return "Normal";
	if (dat == pLogNormal)
		return "LogNormal";
	if (dat == pWeibull)
		return "Weibull";
	if (dat == pTriangle)
		return "Triangle";
	if (dat == pEmpirical)
		return "Empirical";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, prob_distrib_type& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= prob_distrib_type_SIZE)
			THROW("too large of a number\n");
		typeVal = (prob_distrib_type)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < prob_distrib_type_SIZE; ++i)
	{
		typeVal = (prob_distrib_type)i; // casting integer to prob_distrib_type, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading prob_distrib_type\n");
}

//operator for output
ostream& operator<<(ostream& out, prob_distrib_type dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, prob_distrib_type& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

//////
string getName(setStatOp_type dat)
{
	if (dat == sso_none)
		return "none";
	if (dat == sso_valStart)
		return "valStart";
	if (dat == sso_valEnd)
		return "valEnd";
	if (dat == sso_min)
		return "min";
	if (dat == sso_max)
		return "max";
	if (dat == sso_index_min)
		return "index_min";
	if (dat == sso_index_max)
		return "index_max";
	if (dat == sso_mean_arithmetic)
		return "mean_arithmetic";
	if (dat == sso_mean_harmonic)
		return "mean_harmonic";
	if (dat == sso_mean_geometric)
		return "mean_geometric";
	if (dat == sso_sdiv)
		return "sdiv";
	if (dat == sso_var)
		return "var";
	if (dat == sso_cov)
		return "cov";
	if (dat == sso_all_field_val)
		return "all_field_val";
	if (dat == sso_all_field_time_val)
		return "all_field_time_val";
	if (dat == sso_number)
		return "number";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, setStatOp_type& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= setStatOp_type_SIZE)
			THROW("too large of a number\n");
		typeVal = (setStatOp_type)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = -1; i < setStatOp_type_SIZE; ++i)
	{
		typeVal = (setStatOp_type)i; // casting integer to setStatOp_type, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading setStatOp_type\n");
}

//operator for output
ostream& operator<<(ostream& out, setStatOp_type dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, setStatOp_type& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

double getStatvalue(const vector<double>& vals, setStatOp_type sso)
{
	if (sso == sso_cov)
	{
		double mean = getStatvalue(vals, sso_mean_arithmetic);
		double sdiv = getStatvalue(vals, sso_sdiv);
		if (fabs(mean) < 1e-16)
		{
			if (fabs(sdiv) < 1e-16)
				return 0.0;
			return 1e40;
		}
		return sdiv / mean;
	}
	unsigned int sz = vals.size();
	if (sz == 0)
		return 0.0;
	if ((sso == sso_min) || (sso == sso_index_min))
	{
		int index_min = 0;
		double minV = vals[0];
		for (unsigned int i = 0; i < sz; ++i)
			if (vals[i] < minV)
			{
				minV = vals[i];
				index_min = i;
			}
		if (sso == sso_min)
			return minV;
		return (double)index_min;
	}
	if ((sso == sso_max) || (sso == sso_index_max))
	{
		int index_max = 0;
		double maxV = vals[0];
		for (unsigned int i = 0; i < sz; ++i)
			if (vals[i] > maxV)
			{
				maxV = vals[i];
				index_max = i;
			}
		if (sso == sso_max)
			return maxV;
		return (double)index_max;
	}
	if ((sso == sso_mean_arithmetic) || (sso == sso_sdiv) || (sso == sso_var))
	{
		double meanV = 0.0;
		for (unsigned int i = 0; i < sz; ++i)
			meanV += vals[i];
		meanV /= (double)sz;
		if (sso == sso_mean_arithmetic)
			return meanV;
		if (sz <= 1)
			return 0.0;
		double tmp, var = 0.0;

		for (unsigned int i = 0; i < sz; ++i)
		{
			tmp = vals[i] - meanV;
			var += tmp * tmp;
		}
		if (sz > 1)
			var /= (double)(sz - 1);
		if (sso == sso_var)
			return var;
		return sqrt(var);
	}
	if (sso == sso_mean_harmonic)
	{
		double meanV = 0.0;
		for (unsigned int i = 0; i < sz; ++i)
			meanV += 1.0 / vals[i];
		meanV /= (double)sz;
		return 1.0 / meanV;
	}
	if (sso == sso_mean_geometric)
	{
		double meanV = 1.0;
		for (unsigned int i = 0; i < sz; ++i)
			meanV *= vals[i];
		meanV = pow(meanV, 1.0/(double)sz);
		return meanV;
	}
	if (sso == sso_valStart)
		return vals[0];
	if (sso == sso_valEnd)
		return vals[sz - 1];
	if (sso == sso_number)
		return (double)sz;
	cout << sso << '\n';
	THROW("Invalid sso\n");
}

double getWeightedStatvalue(const vector<double>& vals, const vector<double>& weights, setStatOp_type sso)
{
	if (sso == sso_cov)
	{
		double mean = getWeightedStatvalue(vals, weights, sso_mean_arithmetic);
		double sdiv = getWeightedStatvalue(vals, weights, sso_sdiv);
		if (fabs(mean) < 1e-16)
		{
			if (fabs(sdiv) < 1e-16)
				return 0.0;
			return 1e40;
		}
		return sdiv / mean;
	}
	unsigned int sz = vals.size();
	if (sz == 0)
		return 0.0;
	if ((sso == sso_min) || (sso == sso_index_min))
	{
		int index_min = 0;
		double minV = vals[0];
		for (unsigned int i = 0; i < sz; ++i)
			if (vals[i] < minV)
			{
				minV = vals[i];
				index_min = i;
			}
		if (sso == sso_min)
			return minV;
		return (double)index_min;
	}
	if ((sso == sso_max) || (sso == sso_index_max))
	{
		int index_max = 0;
		double maxV = vals[0];
		for (unsigned int i = 0; i < sz; ++i)
			if (vals[i] > maxV)
			{
				maxV = vals[i];
				index_max = i;
			}
		if (sso == sso_max)
			return maxV;
		return (double)index_max;
	}
	if ((sso == sso_mean_arithmetic) || (sso == sso_sdiv) || (sso == sso_var))
	{
		double sumWeight = 0.0, meanV = 0.0;
		for (unsigned int i = 0; i < sz; ++i)
		{
			meanV += weights[i] * vals[i];
			sumWeight += weights[i];
		}
		meanV /= sumWeight;
		if (sso == sso_mean_arithmetic)
			return meanV;
		double tmp, var = 0.0;

		for (unsigned int i = 0; i < sz; ++i)
		{
			tmp = vals[i] - meanV;
			var += weights[i] * tmp * tmp;
		}
		var /= sumWeight;
		if (sso == sso_var)
			return var;
		return sqrt(var);
	}
	if (sso == sso_mean_harmonic)
	{
		double sumWeight = 0.0, meanV = 0.0;
		for (unsigned int i = 0; i < sz; ++i)
		{
			meanV += weights[i] / vals[i];
			sumWeight += weights[i];
		}
		meanV /= sumWeight;
		return 1.0 / meanV;
	}
	if (sso == sso_mean_geometric)
	{
		double sumWeight = 0.0, meanV = 1.0;
		for (unsigned int i = 0; i < sz; ++i)
		{
			meanV *= pow(vals[i], weights[i]);
			sumWeight += weights[i];
		}
		meanV = pow(meanV, 1.0 / sumWeight);
		return meanV;
	}
	if (sso == sso_valStart)
		return vals[0];
	if (sso == sso_number)
		return (double)sz;
	cout << sso << '\n';
	THROW("Invalid sso\n");
}

double getWoWOWeightStatvalue(const vector<double>& vals, bool hasWeight, const vector<double>& weights, setStatOp_type sso)
{
	if (!hasWeight)
		return getStatvalue(vals, sso);
	return getWeightedStatvalue(vals, weights, sso);
}


//////
string getName(valPOperT dat)
{
	if (dat == vpo_none)
		return "none";
	if (dat == vpo_log10)
		return "log10";
	if (dat == vpo_loge)
		return "loge";
	if (dat == vpo_log2)
		return "log2";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string& name, valPOperT& typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= valPOperT_SIZE)
			THROW("too large of a number\n");
		typeVal = (valPOperT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < valPOperT_SIZE; ++i)
	{
		typeVal = (valPOperT)i; // casting integer to valPOperT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading valPOperT\n");
}

//operator for output
ostream& operator<<(ostream& out, valPOperT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

//operator for input
istream& operator>>(istream& in, valPOperT& dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

double get_valPOperT(double valIn, valPOperT oper)
{
	if (oper == vpo_none)
		return valIn;
	if (oper == vpo_log10)
		return log10(valIn);
	if (oper == vpo_loge)
		return log(valIn);
	if (oper == vpo_log2)
		return log2(valIn);
	cout << "oper\t" << oper << '\n';
	THROW("Invalid oper\n");
}

/////

statParamHolder::statParamHolder()
{
	minV = 0.0, maxV = 1.0;
	mode = 0.5;

	mean = 0.0;
	sdiv = 1.0;

	mu =0.0;
	sigma = 1.0;
	
	scale = 1.0;
	shape = 4.0; // m for weibull

	minMax_read = false;
	mode_read = false;
	mean_sdiv_read = false;
	mu_sigma_read = false;
}

void statParamHolder::Read(istream & in) // , double dd2)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "min")
		{
			READ_NDOUBLE(in, buf, minV);
			minMax_read = true;
		}
		else if (buf == "max")
		{
			READ_NDOUBLE(in, buf, maxV);
			minMax_read = true;
		}
		else if (buf == "mode")
		{
			READ_NDOUBLE(in, buf, mode);
			mode_read = true;
		}
		else if (buf == "mean")
		{
			READ_NDOUBLE(in, buf, mean);
			mean_sdiv_read = true;
		}
		else if (buf == "sdiv")
		{
			READ_NDOUBLE(in, buf, sdiv);
			mean_sdiv_read = true;
		}
		else if (buf == "mu")
		{
			READ_NDOUBLE(in, buf, mu);
			mu_sigma_read = true;
		}
		else if (buf == "sigma")
		{
			READ_NDOUBLE(in, buf, sigma);
			mu_sigma_read = true;
		}
		else if (buf == "scale")
		{
			READ_NDOUBLE(in, buf, scale);
		}
		else if (buf == "shape")
		{
			READ_NDOUBLE(in, buf, shape);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
#if 0
	if (dd2 > 0.0)
	{
		double ave = 0.5 * (minV + maxV);
		minV = ave  * (1 - dd2);
		maxV =  ave * (1 + dd2);
		minMax_read = true;
	}
#endif
}

void statParamHolder::WriteData(ostream& out)
{
	out << "{\n";
	if (minMax_read)
	{
		out << "min\n" << minV << '\n';
		out << "max\n" << maxV << '\n';
	}
	if (mode_read)
		out << "mode\n" << mode << '\n';
	if (mean_sdiv_read)
	{
		out << "mean\n" << mean << '\n';
		out << "sdiv\n" << sdiv << '\n';
	}
	if (mu_sigma_read)
	{
		out << "mu\n" << mu << '\n';
		out << "sigma\n" << sigma << '\n';
	}
	out << "scale\n" << scale << '\n';
	out << "shape\n" << shape << '\n';
	out << "}\n";
}

void gRandVar::ReadParameters(istream& in)
{
	paras.Read(in);
	Finalize_Read();
}

double gRandVar::TurnStandardNormalValue2ThisRandom(double standardNormalValue) const
{
	// first calculate standard normal CDF
	double sn_cdf = 0.5 * (1.0 + erf(standardNormalValue / sqrt(2.0)));
	return getInverseCDF(sn_cdf);
}

double gRandVar::getRandomValue() const
{
	double random_01 = slrand() / (double)RAND_MAX;
	// cdf_random_01 = random_01
	return getInverseCDF(random_01);
}

void gRandVar::Finalize_Read()
{
}

gRandVar* gRandVar::CreateCopy()
{
	gRandVar* retVal = NULL;
	if (dist_type == pUniform)
	{
		retVal = new Uniform_RandVar();
		*retVal = *(Uniform_RandVar*)(this);
		return retVal;
	}
	if (dist_type == pTriangle)
	{
		retVal = new Triangle_RandVar();
		*retVal = *(Triangle_RandVar*)(this);
		return retVal;
	}
	if (dist_type == pNormal)
	{
		retVal = new Normal_RandVar();
		*retVal = *(Normal_RandVar*)(this);
		return retVal;
	}
	if (dist_type == pLogNormal)
	{
		retVal = new LogNormal_RandVar();
		*retVal = *(LogNormal_RandVar*)(this);
		return retVal;
	}
	if (dist_type == pWeibull)
	{
		retVal = new Weibull_RandVar();
		*retVal = *(Weibull_RandVar*)(this);
		return retVal;
	}
	if (dist_type == pEmpirical)
	{
		THROW("Not implemented yet\n");
//		retVal = new Empirical_RandVar();
//		*retVal = *(Empirical_RandVar*)(this);
		return retVal;
	}
	cout << dist_type << '\n';
	THROW("copy for this dist_type is not implemented yet\n");
}

void Uniform_RandVar::Finalize_Read()
{
	span = paras.maxV - paras.minV;
	inv_span = 1.0 / span;
}

double Uniform_RandVar::getCDF(double x) const
{
	if (x < paras.minV)
		return 0.0;
	if (x > paras.maxV)
		return 1.0;
	return (x - paras.minV) * inv_span;
}

double Uniform_RandVar::getInverseCDF(double p) const
{
	return paras.minV + span * p;
}

double Uniform_RandVar::getPDF(double x) const
{
	if (x < paras.minV)
		return 0.0;
	if (x > paras.maxV)
		return 0.0;
	return inv_span;
}

double Uniform_RandVar::getRandomValue() const
{
	double random_01 = slrand() / (double)RAND_MAX;
	return paras.minV + span * random_01;
}

double Triangle_RandVar::getCDF(double x) const
{
	if (x < paras.minV)
		return 0.0;
	if (x > paras.maxV)
		return 1.0;
	if (x < paras.mode)
	{
		double tmp = (x - paras.minV);
		return 0.5 * p_md * tmp * tmp * inv_md_m;
	}
	double tmp = (paras.maxV - x);
	return 1.0 - 0.5 * p_md * tmp * tmp * inv_M_md;
}

double Triangle_RandVar::getInverseCDF(double p) const
{
	if (paras.minV - paras.maxV >= -1e-14)
		return paras.mode;
	if (p < cdf_md)
		return sqrt(2.0 * md_m * p / p_md) + paras.minV;
	return paras.maxV - sqrt(2.0 * M_md * (1.0 - p) / p_md);
}

double Triangle_RandVar::getPDF(double x) const
{
	if (x < paras.minV)
		return 0.0;
	if (x > paras.maxV)
		return 0.0;
	if (x < paras.mode)
		return p_md * (x - paras.minV) * inv_md_m;
	return p_md * (paras.maxV - x) * inv_M_md;
}

void Triangle_RandVar::Finalize_Read()
{
	if (!paras.mode_read)
		paras.mode = 0.5 * (paras.minV + paras.maxV);
	span = paras.maxV - paras.minV;
	inv_span = 1.0 / span;
	md_m = paras.mode - paras.minV;
	inv_md_m = 1.0 / md_m;
	M_md = paras.maxV - paras.mode;
	inv_M_md = 1.0 / M_md;
	p_md = 2.0 * inv_span;
	cdf_md = md_m * inv_span;
}

double Normal_RandVar::getCDF(double x) const
{
	return 0.5 * (1.0 + erf((x - paras.mu) / (sqrt(2.0) * paras.sigma)));
}

double Normal_RandVar::getInverseCDF(double p) const
{
	return paras.mean + paras.sdiv *  NormalCDFInverse(p);
}

double Normal_RandVar::getPDF(double x) const
{
	double expv = (x - paras.mu) / paras.sigma;
	return exp(-0.5 * expv * expv) / paras.sigma / sqrt(2.0 * PI);
}

double Normal_RandVar::getRandomValue() const
{
	double r_max = (double)RAND_MAX + 1.0;
	double r1 = (double)slrand();
	double r2 = (double)slrand();
	double v1 = (r1 + 1) / r_max;

	if ((v1 <= 0) || (v1 > 1.0))
		THROW("v1 should between 0 and 1\n");

	return paras.sdiv* sqrt(-2.0 *log(v1))*sin(2.0 * PI * r2 / r_max) + paras.mean;
}

double NormalCDFInverse(double p)
{

	if (p <= 0.0)
	{
		return (double)NUMERIC_NINF;
	}
	if (p >= 1.0)
	{
		return (double)NUMERIC_PINF;
	}

	if (p < 0.5)
	{
		// F^-1(p) = - G^-1(p)
		return -RationalApproximation(sqrt(-2.0*log(p)));
	}
	else
	{
		// F^-1(p) = G^-1(1-p)
		return RationalApproximation(sqrt(-2.0*log(1 - p)));
	}
}

double LogNormal_RandVar::getCDF(double x) const
{
	double log_x = log(x);
	return 0.5 * (1.0 + erf((log_x - paras.mu) / (sqrt(2.0) * paras.sigma)));
}

double LogNormal_RandVar::getInverseCDF(double p) const
{
	double orig_x = paras.mean + paras.sdiv *  NormalCDFInverse(p);
	return exp(orig_x);
}

double LogNormal_RandVar::getPDF(double x) const
{
	double log_x = log(x);
	double expv = (log_x - paras.mu) / paras.sigma;
	return exp(-0.5 * expv * expv) / paras.sigma / sqrt(2.0 * PI);
}

void LogNormal_RandVar::Finalize_Read()
{
	if (paras.mean_sdiv_read)
	{
		double factor = 1 + (paras.sdiv * paras.sdiv / (paras.mean * paras.mean));

		paras.mu = log(paras.mean / sqrt(factor));
		paras.sigma = sqrt(log(factor));
	}
}

double Weibull_RandVar::getCDF(double x) const
{
	if (x < paras.minV)
		return 0.0;
	double del_x = x - paras.minV;
	double expv = del_x / paras.scale;
	expv = pow(expv, paras.shape);
	return 1.0 - exp(-expv);
}

double Weibull_RandVar::getInverseCDF(double p) const
{
	if (1.0 - p < 1e-7)
		return NUMERIC_PINF;
	double alpha = -log(1.0 - p); // alpha = (x/scale)^shape
	return paras.minV + paras.scale * pow(alpha, 1.0 / paras.shape);
}

double Weibull_RandVar::getPDF(double x) const
{
	if (x < paras.minV)
		return 0.0;
	double del_x = x - paras.minV;
	double scaledVal = del_x / paras.scale;
	double expv = pow(scaledVal, paras.shape - 1.0);
	return paras.shape / paras.scale * expv * exp(-(scaledVal * expv));
}

double Weibull_RandVar::getRandomValue() const
{
	double r = 1.0;

	double r_max = (double)RAND_MAX + 1.0;
	double r1 = (double)slrand();
	double r2 = (double)slrand();
	double y = (r1 + 1) / r_max;
	static double tol = 1e-10;

	if ((y <= 0) || (y > 1.0))
		THROW("y should between 0 and 1\n");

	double omy = 1.0 - y;
	if (omy < tol)
		omy = tol;
	return paras.minV + pow(-1.0 / r * log(omy), 1.0 / paras.shape) * paras.scale;
}

gRandVar * gRandVar_Factory(prob_distrib_type dist_type)
{
	gRandVar* retVal = NULL;
	if (dist_type == pUniform)
	{
		retVal = new Uniform_RandVar();
	}
	else if (dist_type == pTriangle)
	{
		retVal = new Triangle_RandVar();
	}
	else if (dist_type == pNormal)
	{
		retVal = new Normal_RandVar();
	}
	else if (dist_type == pLogNormal)
	{
		retVal = new LogNormal_RandVar();
	}
	else if (dist_type == pWeibull)
	{
		retVal = new Weibull_RandVar();
	}
	else if (dist_type == pEmpirical)
	{
		THROW("Function option is not implemented yet\n");
	}
	if (retVal != NULL)
		retVal->dist_type = dist_type;
	return retVal;
}

gRandVar* ReadRandomVariableFromFile(istream& in)
{
	prob_distrib_type dist_type;
	in >> dist_type;
	gRandVar * retVal = gRandVar_Factory(dist_type);
	if (retVal != NULL)
		retVal->ReadParameters(in);
	return retVal;
}

double RationalApproximation(double t)
{
	// Abramowitz and Stegun formula 26.2.23.
	// The absolute value of the error should be less than 4.5 e-4.
	double c[] = { 2.515517, 0.802853, 0.010328 };
	double d[] = { 1.432788, 0.189269, 0.001308 };
	return t - ((c[2] * t + c[1])*t + c[0]) /
		(((d[2] * t + d[1])*t + d[0])*t + 1.0);
}

inline int slrand()
{
#if 0
	if (rand_cntr0 < ULONG_MAX)
		++rand_cntr0;
	else
	{
		rand_cntr0 = 0;
		++rand_cntr1;
	}
#endif
	return rand();
}
