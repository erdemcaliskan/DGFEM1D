#ifndef COMMON_MACROS__H
#define COMMON_MACROS__H

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#include "FileIOMacros.h"

#ifndef TEMPORARY_DISABLE_SNIPPET
#define TEMPORARY_DISABLE_SNIPPET 1
#endif

#ifndef TENSOR_KEY_OLD_VERSION
#define TENSOR_KEY_OLD_VERSION 0
#endif

// these flags force the code IO (for simplex geometry + polynomial data / other
// futher data to be ASCII binary)
#define FORCE_CIO_ASCII 0
#define FORCE_CIO_BINARY 0

// in the following isBinary is ignored is either of the FORCE flags above are 1

#if FORCE_CIO_BINARY
#define READ_CIO_INTEGER(in, dat, isBinary) in.read((char *)&dat, sizeof(int));
#else
#if FORCE_CIO_ASCII
#define READ_CIO_INTEGER(in, dat, isBinary) in >> dat;
#else
#define READ_CIO_INTEGER(in, dat, isBinary)                                                                            \
    {                                                                                                                  \
        if (isBinary)                                                                                                  \
            in.read((char *)&dat, sizeof(int));                                                                        \
        else                                                                                                           \
            in >> dat;                                                                                                 \
    }
#endif
#endif

#if FORCE_CIO_BINARY
#define READ_CIO_DOUBLE(in, dat, isBinary) in.read((char *)&dat, sizeof(double));
#else
#if FORCE_CIO_ASCII
#define READ_CIO_DOUBLE(in, dat, isBinary) in >> dat;
#else
#define READ_CIO_DOUBLE(in, dat, isBinary)                                                                             \
    {                                                                                                                  \
        if (isBinary)                                                                                                  \
            in.read((char *)&dat, sizeof(double));                                                                     \
        else                                                                                                           \
            in >> dat;                                                                                                 \
    }
#endif
#endif

#if FORCE_CIO_BINARY
#define WRITE_CIO_INTEGER(out, dat, isBinary) out.write((char *)&dat, sizeof(int));
#else
#if FORCE_CIO_ASCII
#define WRITE_CIO_INTEGER(out, dat, isBinary) out << dat;
#else
#define WRITE_CIO_INTEGER(out, dat, isBinary)                                                                          \
    {                                                                                                                  \
        if (isBinary)                                                                                                  \
            out.write((char *)&dat, sizeof(int));                                                                      \
        else                                                                                                           \
            out << dat;                                                                                                \
    }
#endif
#endif

#if FORCE_CIO_BINARY
#define WRITE_CIO_DOUBLE(out, dat, isBinary) out.write((char *)&dat, sizeof(double));
#else
#if FORCE_CIO_ASCII
#define WRITE_CIO_DOUBLE(out, dat, isBinary) out << dat;
#else
#define WRITE_CIO_DOUBLE(out, dat, isBinary)                                                                           \
    {                                                                                                                  \
        if (isBinary)                                                                                                  \
            out.write((char *)&dat, sizeof(double));                                                                   \
        else                                                                                                           \
            out << dat;                                                                                                \
    }
#endif
#endif

#if FORCE_CIO_BINARY
#define WRITE_CIO_TAB(out, isBinary)
#else
#if FORCE_CIO_ASCII
#define WRITE_CIO_TAB(out, isBinary) out << '\t';
#else
#define WRITE_CIO_TAB(out, isBinary)                                                                                   \
    {                                                                                                                  \
        if (!isBinary)                                                                                                 \
            out << '\t';                                                                                               \
    }
#endif
#endif

#if FORCE_CIO_BINARY
#define WRITE_CIO_EOL(out, isBinary)
#else
#if FORCE_CIO_ASCII
#define WRITE_CIO_EOL(out, isBinary) out << '\n';
#else
#define WRITE_CIO_EOL(out, isBinary)                                                                                   \
    {                                                                                                                  \
        if (!isBinary)                                                                                                 \
            out << '\n';                                                                                               \
    }
#endif
#endif

template <class T> inline void write_cio_midline(std::ostream &out, const T &dat, bool isBinary = true)
{
#if FORCE_CIO_BINARY
    out.write((char *)&dat, sizeof(T));
    return;
#endif
#if FORCE_CIO_ASCII
    out << dat << '\t';
    return;
#endif
    if (isBinary)
    {
        out.write((char *)&dat, sizeof(T));
        return;
    }
    else
    {
        out << dat << '\t';
        return;
    }
}

template <class T> inline void write_cio_endline(std::ostream &out, const T &dat, bool isBinary = true)
{
#if FORCE_CIO_BINARY
    out.write((char *)&dat, sizeof(T));
    return;
#endif
#if FORCE_CIO_ASCII
    out << dat << '\n';
    return;
#endif
    if (isBinary)
    {
        out.write((char *)&dat, sizeof(T));
        return;
    }
    else
    {
        out << dat << '\n';
        return;
    }
}

#define to_lower(str)                                                                                                  \
    {                                                                                                                  \
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);                                                \
    }

#define LOWER_CASE(str)                                                                                                \
    {                                                                                                                  \
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);                                                \
    }

template <class T> bool fromString(const string &name, T &dat)
{
    istringstream ss(name);
    ss >> dat;
    if (ss.fail())
        return false;
    return true;
}

template <class T> void toString(T &dat, string &name)
{
    ostringstream convert;
    convert << dat;
    name = convert.str();
}

// template <class Ty> istream &binary_read(istream &in, string &buf, Ty &var) {
//   if (!in.read((char *)&var, sizeof(Ty))) {
//     stringstream ss;
//     ss << "Error occured during binary file read::number of bytes read["
//        << in.gcount() << "]\n";
//     THROW(ss.str().c_str());
//   }
//   buf = to_string(var);
//   return in;
// }
//
// template <class Ty> ostream &binary_write(ostream &out, Ty var) {
//   if (!out.write((char *)&var, sizeof(Ty))) {
//     THROW("Error occured during binary file write");
//   }
//   return out;
// }

// reads the next input, if it's not the write type, it goes back and returns
// false if the next input starts with character # the entire line is read

// readUntilSuccessful
//	true	reads until it finds specific type requested and skips all
// incorrect formats and comments 			only returns false if
// until the end of the file the particular type cannot be found 	false
// reads ONLY one value: 				if it has correct format
// =>		returns true if NOT
//=>		puts back the wrong value into string
// comment ends with precedingEndCommend commentChar, for example for
// commentChar = '#' and precedingEndCommend = '\', end of comment would be \#
// 0 means unacceptable format
// -1 means end of file (fail to read)
// 1 means successful read
template <class T>
inline int ReadSpecificType(istream &in, T &val, string &buf, bool readUntilSuccessful = false, char commentChar = '#',
                            string precedingEndCommend = "\\#")
{
    string blockLeftBrac = "%{";
    string blockRightBrac = "%}";

    in >> buf;
    if (in.fail() == true)
    {
        buf = "";
        return -1;
    }
    size_t sz = buf.length();
    if ((sz > 0) && (buf[0] == commentChar))
    {
        char c;
        in.get(c);
        bool brk = false;
        while (c != '\n')
        {
            in.get(c);
            if (c == precedingEndCommend[0])
            {
                in.putback(c);
                in >> buf;
                if (buf == precedingEndCommend)
                    brk = true;
            }
            if (brk)
                break;
        }
        //		getline(in, buf);
        return ReadSpecificType(in, val, buf, readUntilSuccessful, commentChar);
    }
    // BLOCK READ -- PLC
    if ((sz >= blockLeftBrac.length()) && buf.compare(0, blockLeftBrac.length(), blockLeftBrac) == 0)
    {
        char c;
        in.get(c);
        bool brk = false;
        while (!brk)
        {
            in.get(c);
            if (c == blockRightBrac[0])
            {
                in.putback(c);
                in >> buf;
                if (buf == blockRightBrac)
                    brk = true;
            }
        }
        //		getline(in, buf);
        return ReadSpecificType(in, val, buf, readUntilSuccessful, commentChar);
    }
    if (fromString(buf, val) == false)
    {
        if (readUntilSuccessful == true)
            return ReadSpecificType(in, val, buf, readUntilSuccessful, commentChar);
        for (size_t i = 0; i < sz; ++i)
            in.putback(buf[sz - i - 1]);
        return 0;
    }
    return 1;
}

#define READ_INTEGER(iFile, buf, var)                                                                                  \
    {                                                                                                                  \
        string bufs;                                                                                                   \
        if (ReadSpecificType(iFile, var, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading integer\n");                                                        \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_OINTEGER(iFile, buf, var)                                                                                 \
    {                                                                                                                  \
        string bufs;                                                                                                   \
        if (ReadSpecificType(iFile, var, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading integer\n");                                                        \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_NINTEGER(iFile, bufs, var)                                                                                \
    {                                                                                                                  \
        if (ReadSpecificType(iFile, var, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << bufs << '\n';                                                                           \
            THROW("invalid file format for reading integer\n");                                                        \
        }                                                                                                              \
    }

#define READ_BOOL(iFile, buf, var)                                                                                     \
    {                                                                                                                  \
        string tmp;                                                                                                    \
        string bufs;                                                                                                   \
        if (ReadSpecificType(iFile, tmp, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading boolean\n");                                                        \
        }                                                                                                              \
        if ((tmp[0] == 'y') || (tmp[0] == '1'))                                                                        \
            var = true;                                                                                                \
        else if ((tmp[0] == 'n') || (tmp[0] == '0'))                                                                   \
            var = false;                                                                                               \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_OBOOL(iFile, buf, var)                                                                                    \
    {                                                                                                                  \
        string tmp;                                                                                                    \
        string bufs;                                                                                                   \
        if (ReadSpecificType(iFile, tmp, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading boolean\n");                                                        \
        }                                                                                                              \
        if ((tmp[0] == 'y') || (tmp[0] == '1'))                                                                        \
            var = true;                                                                                                \
        else if ((tmp[0] == 'n') || (tmp[0] == '0'))                                                                   \
            var = false;                                                                                               \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_NBOOL(iFile, bufs, var)                                                                                   \
    {                                                                                                                  \
        string tmp;                                                                                                    \
        if (ReadSpecificType(iFile, tmp, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << bufs << '\n';                                                                           \
            THROW("invalid file format for reading boolean\n");                                                        \
        }                                                                                                              \
        if ((tmp[0] == 'y') || (tmp[0] == '1'))                                                                        \
            var = true;                                                                                                \
        else if ((tmp[0] == 'n') || (tmp[0] == '0'))                                                                   \
            var = false;                                                                                               \
    }

#define READ_ODOUBLE(iFile, buf, var)                                                                                  \
    {                                                                                                                  \
        string bufs;                                                                                                   \
        if (ReadSpecificType<double>(iFile, var, bufs) == 0)                                                           \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading double\n");                                                         \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_NDOUBLE(iFile, bufs, var)                                                                                 \
    {                                                                                                                  \
        if (ReadSpecificType<double>(iFile, var, bufs) == 0)                                                           \
        {                                                                                                              \
            cout << "buf\n" << bufs << '\n';                                                                           \
            THROW("invalid file format for reading double\n");                                                         \
        }                                                                                                              \
    }

#define READ_NSTRING(iFile, buf, var)                                                                                  \
    {                                                                                                                  \
        if (ReadSpecificType(iFile, var, buf) == 0)                                                                    \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading double\n");                                                         \
        }                                                                                                              \
    }

#define READ_OSTRING(iFile, buf, var)                                                                                  \
    {                                                                                                                  \
        string bufs, s;                                                                                                \
        if (ReadSpecificType(iFile, s, bufs) == 0)                                                                     \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading double\n");                                                         \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
        strcpy(var, s.c_str());                                                                                        \
    }

#endif

#define READ_OCHAR(iFile, buf, var)                                                                                    \
    {                                                                                                                  \
        string bufs;                                                                                                   \
        if (ReadSpecificType(iFile, var, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading char\n");                                                           \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_NCHAR(iFile, bufs, var)                                                                                   \
    {                                                                                                                  \
        if (ReadSpecificType(iFile, var, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << bufs << '\n';                                                                           \
            THROW("invalid file format for reading char\n");                                                           \
        }                                                                                                              \
    }

#if CODE_FLAG != OLD_CODE

#define READ_BOOL(iFile, buf, var)                                                                                     \
    {                                                                                                                  \
        string tmp;                                                                                                    \
        string bufs;                                                                                                   \
        if (ReadSpecificType(iFile, tmp, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading boolean\n");                                                        \
        }                                                                                                              \
        if ((tmp[0] == 'y') || (tmp[0] == '1'))                                                                        \
            var = true;                                                                                                \
        else if ((tmp[0] == 'n') || (tmp[0] == '0'))                                                                   \
            var = false;                                                                                               \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_DOUBLE(iFile, buf, var)                                                                                   \
    {                                                                                                                  \
        string bufs;                                                                                                   \
        if (ReadSpecificType<double>(iFile, var, bufs) == 0)                                                           \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading double\n");                                                         \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#define READ_STRING(iFile, buf, var)                                                                                   \
    {                                                                                                                  \
        string bufs, s;                                                                                                \
        if (ReadSpecificType(iFile, s, bufs) == 0)                                                                     \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading double\n");                                                         \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
        strcpy(var, s.c_str());                                                                                        \
    }

#define READ_CHAR(iFile, buf, var)                                                                                     \
    {                                                                                                                  \
        string bufs;                                                                                                   \
        if (ReadSpecificType(iFile, var, bufs) == 0)                                                                   \
        {                                                                                                              \
            cout << "buf\n" << buf << '\n';                                                                            \
            THROW("invalid file format for reading char\n");                                                           \
        }                                                                                                              \
        strcpy(buf, bufs.c_str());                                                                                     \
    }

#endif

//------------------------------------------------------------------------------------
// intermediate implementation (until 11/2014)
#if 0
#define READ_INTEGER(iFile, buf, var)                                                                                  \
    {                                                                                                                  \
        while (iFile.getline(buf, sizeof(buf)) != NULL)                                                                \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == '\n') || (buf[0] == '\0'))                                               \
                continue;                                                                                              \
            else                                                                                                       \
            {                                                                                                          \
                if (sscanf(buf, "%d", &var) <= 0)                                                                      \
                    continue;                                                                                          \
                else                                                                                                   \
                    break;                                                                                             \
            }                                                                                                          \
        }                                                                                                              \
    }

#define READ_BOOL(iFile, buf, var)                                                                                     \
    {                                                                                                                  \
        while (iFile.getline(buf, sizeof(buf)) != NULL)                                                                \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == '\n') || (buf[0] == '\0'))                                               \
                continue;                                                                                              \
            else                                                                                                       \
            {                                                                                                          \
                if (sscanf(buf, "%s", &buf) <= 0)                                                                      \
                    continue;                                                                                          \
                else                                                                                                   \
                    break;                                                                                             \
            }                                                                                                          \
        }                                                                                                              \
        if ((buf[0] == 'y') || (buf[0] == '1'))                                                                        \
            var = true;                                                                                                \
        else if ((buf[0] == 'n') || (buf[0] == '0'))                                                                   \
            var = false;                                                                                               \
        else                                                                                                           \
        {                                                                                                              \
            std::cout << "invalid boolean input!\n";                                                                   \
            getchar();                                                                                                 \
            getchar();                                                                                                 \
            exit(0);                                                                                                   \
        }                                                                                                              \
    }

#define READ_DOUBLE(iFile, buf, var)                                                                                   \
    {                                                                                                                  \
        while (iFile.getline(buf, sizeof(buf)) != NULL)                                                                \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == '\n') || (buf[0] == '\0'))                                               \
                continue;                                                                                              \
            else                                                                                                       \
            {                                                                                                          \
                if (sscanf(buf, "%s", &buf) <= 0)                                                                      \
                    continue;                                                                                          \
                else                                                                                                   \
                    break;                                                                                             \
            }                                                                                                          \
        }                                                                                                              \
        var = atof(buf);                                                                                               \
    }

#define READ_STRING(iFile, buf, var)                                                                                   \
    {                                                                                                                  \
        while (iFile.getline(buf, sizeof(buf)) != NULL)                                                                \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == '\n') || (buf[0] == '\0'))                                               \
                continue;                                                                                              \
            else                                                                                                       \
            {                                                                                                          \
                if (sscanf(buf, "%s", &var) <= 0)                                                                      \
                    continue;                                                                                          \
                else                                                                                                   \
                    break;                                                                                             \
            }                                                                                                          \
        }                                                                                                              \
    }

#define READ_INTEGER_ARRAY(iFile, buf, var)                                                                            \
    {                                                                                                                  \
        while (iFile.getline(buf, sizeof(buf)) != NULL)                                                                \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == '\n') || (buf[0] == '\0'))                                               \
                continue;                                                                                              \
            else                                                                                                       \
                break;                                                                                                 \
        }                                                                                                              \
        char *par;                                                                                                     \
        int nRows = strtol(buf, &par, 10);                                                                             \
        if (!((*par == ' ') || (*par == '\n') || (*par == '\0')))                                                      \
        {                                                                                                              \
            cout << "Error: Should be integer in integer_array:" << endl;                                              \
            exit(1);                                                                                                   \
        }                                                                                                              \
        for (int j = 0; j < nRows; j++)                                                                                \
        {                                                                                                              \
            iFile.getline(buf, sizeof(buf));                                                                           \
            char *pch = strtok(buf, " ");                                                                              \
            int k = 0;                                                                                                 \
            while ((pch != NULL))                                                                                      \
            {                                                                                                          \
                var[j][k] = strtol(pch, &par, 10);                                                                     \
                k++;                                                                                                   \
                if (!((*par == ' ') || (*par == '\n') || (*par == '\0')))                                              \
                {                                                                                                      \
                    cout << "Error: Should be Integer in integer_array:" << endl;                                      \
                    exit(1);                                                                                           \
                }                                                                                                      \
                pch = strtok(NULL, " ");                                                                               \
            }                                                                                                          \
        }                                                                                                              \
    }

#define READ_DOUBLE_ARRAY(iFile, buf, var)                                                                             \
    {                                                                                                                  \
        while (iFile.getline(buf, sizeof(buf)) != NULL)                                                                \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == '\n') || (buf[0] == '\0'))                                               \
                continue;                                                                                              \
            else                                                                                                       \
                break;                                                                                                 \
        }                                                                                                              \
        char *pad;                                                                                                     \
        int dRows = strtol(buf, &pad, 10);                                                                             \
        if (!((*pad == ' ') || (*pad == '\n') || (*pad == '\0')))                                                      \
        {                                                                                                              \
            cout << "Error: Should be integer in double_array:" << endl;                                               \
            exit(1);                                                                                                   \
        }                                                                                                              \
        for (int j = 0; j < dRows; j++)                                                                                \
        {                                                                                                              \
            iFile.getline(buf, sizeof(buf));                                                                           \
            char *pdh = strtok(buf, " ");                                                                              \
            int k = 0;                                                                                                 \
            while ((pdh != NULL))                                                                                      \
            {                                                                                                          \
                var[j][k] = strtod(pdh, &pad);                                                                         \
                k++;                                                                                                   \
                if (!((*pad == ' ') || (*pad == '\n') || (*pad == '\0')))                                              \
                {                                                                                                      \
                    cout << "Error: Should be double in double_array:" << endl;                                        \
                    exit(1);                                                                                           \
                }                                                                                                      \
                pdh = strtok(NULL, " ");                                                                               \
            }                                                                                                          \
        }                                                                                                              \
    }

#define READ_CHAR(iFile, buf, var)                                                                                     \
    {                                                                                                                  \
        while (iFile.getline(buf, sizeof(buf)) != NULL)                                                                \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == '\n') || (buf[0] == '\0'))                                               \
                continue;                                                                                              \
            else                                                                                                       \
                break;                                                                                                 \
        }                                                                                                              \
        var = buf[0];                                                                                                  \
    }
#endif

// old ones
#if 0 //	old old implementation
#define READ_INT(fp, buf, var)                                                                                         \
    {                                                                                                                  \
        while (fgets(buf, sizeof(buf), fp) != NULL)                                                                    \
        {                                                                                                              \
            if ((buf[0] == '#') || (buf[0] == ' ') || (buf[0] == '\n'))                                                \
                continue;                                                                                              \
            else                                                                                                       \
                break;                                                                                                 \
        }                                                                                                              \
        var = atoi(buf);                                                                                               \
    }

#define FPREAD_STRING(fp, buf, var)                                                                                    \
    while (fgets(buf, sizeof(buf), fp) != NULL)                                                                        \
    {                                                                                                                  \
        if ((buf[0] == '#') || (buf[0] == ' ') || (buf[0] == '\n'))                                                    \
            continue;                                                                                                  \
        else                                                                                                           \
            break;                                                                                                     \
    }                                                                                                                  \
    sscanf(buf, "%s", var);
#define FPREAD_INT(fp, buf, var)                                                                                       \
    while (fgets(buf, sizeof(buf), fp) != NULL)                                                                        \
    {                                                                                                                  \
        if ((buf[0] == '#') || (buf[0] == ' ') || (buf[0] == '\n'))                                                    \
            continue;                                                                                                  \
        else                                                                                                           \
            break;                                                                                                     \
    }                                                                                                                  \
    var = atoi(buf);
#define FPREAD_DOUBLE(fp, buf, var)                                                                                    \
    while (fgets(buf, sizeof(buf), fp) != NULL)                                                                        \
    {                                                                                                                  \
        if ((buf[0] == '#') || (buf[0] == ' ') || (buf[0] == '\n'))                                                    \
            continue;                                                                                                  \
        else                                                                                                           \
            break;                                                                                                     \
    }                                                                                                                  \
    var = atof(buf);
#define FPREAD_CHAR(fp, buf, var)                                                                                      \
    while (fgets(buf, sizeof(buf), fp) != NULL)                                                                        \
    {                                                                                                                  \
        if ((buf[0] == '#') || (buf[0] == ' ') || (buf[0] == '\n'))                                                    \
            continue;                                                                                                  \
        else                                                                                                           \
            break;                                                                                                     \
    }                                                                                                                  \
    var = buf[0];                                                                                                      \
#endif

#endif