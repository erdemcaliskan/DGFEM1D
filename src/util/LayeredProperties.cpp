#include "LayeredProperties.h"
#include "commonMacros.h"

istream & operator>>(istream & in, oneBulk_Elastic_Prop & dat)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return in;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if ((buf == "E") || (buf == "E_iso"))
		{
			READ_NINTEGER(in, buf, dat.E);
		}
		else if ((buf == "rho") || (buf == "r"))
		{
			READ_NINTEGER(in, buf, dat.rho);
		}
		else if ((buf == "d") || (buf == "damping") || (buf == "D_vv"))
		{
			READ_NINTEGER(in, buf, dat.damping);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	return in;
}


void oneBulk_Elastic_Prop::setFinalValue(const oneBulk_Elastic_Prop& baseMat, const Bulk_Elastic_Modifier* bemPtr)
{
	E = baseMat.E * bemPtr->CFactor;
	rho = baseMat.rho * bemPtr->rhoFactor;
	damping = baseMat.damping * bemPtr->dampingFactor;
	length = bemPtr->length;
}

oneBulk_Elastic_Prop::oneBulk_Elastic_Prop()
{
	E = 1.0;
	rho = 1.0;
	damping = 0.0;
}

Bulk_Elastic_Modifier::Bulk_Elastic_Modifier()
{
	length = 1.0;
	bulk_flag = 1;
	CFactor = 1.0, rhoFactor = 1.0, dampingFactor = 1.0;
	b_modifies = false;
}

bool Bulk_Elastic_Modifier::Read_Bulk_Elastic_Modifier(istream& in)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
	{
		if (buf == "}")
			return false;
		else
		{
			THROW("istream should start with {");
		}
	}
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		if (buf == "l")
		{
			READ_NDOUBLE(in, buf, length);
		}
		else if (buf == "f")
		{
			READ_NINTEGER(in, buf, bulk_flag);
		}
		else if (buf == "C")
		{
			READ_NDOUBLE(in, buf, CFactor);
		}
		else if (buf == "r")
		{
			READ_NDOUBLE(in, buf, rhoFactor);
		}
		else if (buf == "d")
		{
			READ_NDOUBLE(in, buf, dampingFactor);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	Initialize_Bulk_Elastic_Modifier();
	return true;
}

void Bulk_Elastic_Modifier::Print(ostream & out) const
{
	out << "bulk_flag\t" << bulk_flag << '\t';
	out << "length\t" << length << '\t';
	out << "CFactor\t" << CFactor << '\t';
	out << "rhoFactor\t" << rhoFactor << '\t';
	out << "dampingFactor\t" << dampingFactor << '\t';
	out << "b_modifies\t" << b_modifies << '\n';
}

void Bulk_Elastic_Modifier::Initialize_Bulk_Elastic_Modifier()
{
	b_modifies = false;
	static double tol = 1e-7;
	if (fabs(CFactor - 1.0) > tol)
	{
		b_modifies = true;
		return;
	}
	if (fabs(rhoFactor - 1.0) > tol)
	{
		b_modifies = true;
		return;
	}
	if (fabs(dampingFactor - 1.0) > tol)
	{
		dampingFactor = true;
		b_modifies = true;
		return;
	}
}

Bulk_Elastic_Prop::Bulk_Elastic_Prop()
{
	directSpaceSizeModifier = 0.0;
	numRepeatSequence = 1;
	b_directSpaceSizeModifier = false;
}

void Bulk_Elastic_Prop::Read_Bulk_Elastic_Prop(istream& in, int serialNumber, bool readData)
{
	string serialNumber_str = "";
	if (serialNumber >= 0)
	{
		toString(serialNumber, serialNumber_str);
		serialNumber_str = "_" + serialNumber_str;
	}
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
		if (buf == "numRepeatSequence")
		{
			READ_NINTEGER(in, buf, numRepeatSequence);
		}
		else if (buf == "directSpaceSizeModifier")
		{
			READ_NDOUBLE(in, buf, directSpaceSizeModifier);
			b_directSpaceSizeModifier = (fabs(directSpaceSizeModifier) > DBL_MIN);
		}
		else if (buf == "bulk_elastic_map")
		{
			READ_NSTRING(in, buf, buf);
			if (buf != "{")
			{
				THROW("start of block should start with {");
			}
			READ_NSTRING(in, buf, buf);
			while (buf != "}")
			{
				GID id;
				if (fromString(buf, id) == false)
				{
					cout << "buf\t" << buf << "invalid format\n";
				}
				in >> baseBulkProperties[id];
				READ_NSTRING(in, buf, buf);
			}
		}
		else if (buf == "bulkModifierFileNameWOserExt")
		{
			string fileName4Modifier;
			READ_NSTRING(in, buf, fileName4Modifier);
			if (readData)
			{
				fstream indat;
#if 0
				if (fileName4Modifier == "inline")
				{
					indat.open(fileName.c_str(), ios::in);
					buf = "none";
					while (buf != "bulkModifiers_direct")
						READ_NSTRING(in, buf, buf);
				}
				else
#endif
				{
					fileName4Modifier += serialNumber_str;
					fileName4Modifier += ".txt";
					indat.open(fileName4Modifier.c_str(), ios::in);
					if (!indat.is_open())
					{
						cout << "fileName4Modifier\t" << fileName4Modifier << '\n';
						THROW("Cannot open file\n");
					}
				}
				READ_NSTRING(indat, buf, buf);
				if (buf != "{")
				{
					THROW("start of block should start with {");
				}
				Bulk_Elastic_Modifier bem;
				bool cont_reading = bem.Read_Bulk_Elastic_Modifier(indat);
				while (cont_reading)
				{
					bulkModifiers.push_back(bem);
					cont_reading = bem.Read_Bulk_Elastic_Modifier(indat);
				}
			}
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
	InitializeAfterRead();
}

void Bulk_Elastic_Prop::InitializeAfterRead()
{
	bool min_bulk_size = (directSpaceSizeModifier > 1e-40);
	unsigned int num_refinement = 1;
	bool b_num_refinement = (directSpaceSizeModifier < -1);
	bool checkSize = (min_bulk_size || b_num_refinement);
	if (b_num_refinement)
		num_refinement = MAX((int)round(-directSpaceSizeModifier), 1);

	finalBulkSegments.clear();
	unsigned int num_modifiers = bulkModifiers.size();
	double current_length, new_length;
	for (unsigned int i = 0; i < num_modifiers; ++i)
	{
		Bulk_Elastic_Modifier bem = bulkModifiers[i];
		current_length = bem.length;
		if (b_directSpaceSizeModifier)
		{
			if (min_bulk_size)
				num_refinement = (int)ceil(current_length / directSpaceSizeModifier);
		}
		if (num_refinement > 1)
		{
			new_length = current_length / num_refinement;
			Bulk_Elastic_Modifier bemNew = bem;
			bemNew.length = new_length;
			oneBulk_Elastic_Prop bep;
			bep.setFinalValue(baseBulkProperties[bem.bulk_flag], &bemNew);
			for (unsigned int i = 0; i < num_refinement; ++i)
				finalBulkSegments.push_back(bep);
		}
		else
		{
			oneBulk_Elastic_Prop bep;
			bep.setFinalValue(baseBulkProperties[bem.bulk_flag], &bem);
			finalBulkSegments.push_back(bep);
		}
	}
	sz_oneSequence = finalBulkSegments.size();
	sz_allsequences = sz_oneSequence * numRepeatSequence;
	if (numRepeatSequence == 1)
		return;
	finalBulkSegments.resize(sz_allsequences);
	unsigned int st = sz_oneSequence;
	for (unsigned int si = 1; si < numRepeatSequence; ++si)
	{
		for (unsigned int i = 0; i < sz_oneSequence; ++i)
			finalBulkSegments[st++] = finalBulkSegments[i];
	}
}
