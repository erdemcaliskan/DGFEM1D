/**
FILEIOMacros.h
Purpose:  contains the macro definitions for file management that is standard
across OS platforms.

@version 1.0  6/13/2016
*/

#ifndef FILE_IO_MACROS__H
#define FILE_IO_MACROS__H

/** C++ Standard and STL includes*/
//=========================
#include <cstddef>
#include <cstring> //hang revised 2018.1.22 #include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
//#include <cstdlib>
#include <stdio.h>
//#include <cstdio>
#include <cerrno>
#include <ctime>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef VCPP
#define BKSLSH "\\"
#else
#define BKSLSH "/"
#endif

#define SUPP_OUT 1

#if SUPP_OUT
#ifndef cmdSuffix
#define cmdSuffix " >nul 2>nul"
#endif
#else
#ifndef cmdSuffix
#define cmdSuffix ""
#endif
#endif

using namespace std;

#ifdef _MSC_VER
#include <direct.h>
#include <windows.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#ifndef VCPP
template <class Ty> string to_string(Ty val) {
  string strVal;
  ostringstream oss;
  oss << val;
  strVal = oss.str();
  return strVal;
}
#endif

/*IO FUNCTIONS*/
// Peek stream for next character
inline char PEEK(istream &in) {
  string buf;
  char out, temp;
  in.get(temp);

  while (isspace(temp) || temp == '\n' || temp == '#') {
    if (temp == '#')
      getline(in, buf);

    in.get(temp);
  }

  in.putback(temp);
  out = temp;
  return out;
}

template <class T>
int vectorOut(ostream &out, const vector<T> &dat, string delim = "{}") {
  out << delim.at(0) << " ";
  for (int i = 0; i < (int)dat.size(); ++i)
    out << dat[i] << " ";
  out << delim.at(1);

  return (int)dat.size();
}

typedef enum {
  appendF,
  renameF,
  copyF,
  moveF,
  removeD,
  makeD,
  removeF
} op; // F:File, D:Directory

inline bool CMD(const char *command) // command line operation
{
  int stat;
  string cmd;
#if defined(_MSC_VER) // implemented in case future implementation can be
                      // optimized
  cmd = "\"" + (string)command + cmdSuffix + "\"";
  stat = system(cmd.c_str());
#else
  cmd = (string)command + cmdSuffix;
  stat = system(cmd.c_str());
#endif

  return (bool)stat;
}

inline string getFileNameFromPath(string str) {
  size_t found = str.find_last_of("//\\/");
  string filename = str.substr(found + 1);
  if (filename.empty())
    filename = str;
  return filename;
}

inline vector<string> getFileParts(const string &str) {
  vector<string> parts(3);
  size_t found = str.find_last_of("//\\/");
  string filePath = str.substr(0, found);
  string filename = str.substr(found + 1);
  string extension = "";
  if (!filename.empty()) {
    found = filename.find_last_of(".");
    extension = filename.substr(found + 1);
    filename = filename.substr(0, found);
  }
  parts[0] = filePath;
  parts[1] = filename;
  parts[2] = extension;

  return parts;
}

inline int getCurrWorkDir(string &str) {
  char cCurrentPath[FILENAME_MAX];

  if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) {
    return errno;
  }

  cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

  str = (string)cCurrentPath;
  return 1;
}

inline int path_attribute(const char *path) {
  string prePath;
  getCurrWorkDir(prePath);
  string tempPath = string(path);

  string newPath;
  size_t found = tempPath.find(prePath);
  if (found != std::string::npos) {
    newPath = tempPath;
  } else {
    newPath = prePath + BKSLSH + tempPath;
  }

  found = newPath.find_last_of("//\\/");
  string filename = newPath.substr(found + 1);
  if (filename.empty())
    newPath = newPath.substr(0, found);

  struct stat
      s_buf; // ref:
             // http://pubs.opengroup.org/onlinepubs/007908799/xsh/sysstat.h.html
  int status = stat(newPath.c_str(), &s_buf);

  if (status == 0) // the path exist and is something
  {
#if defined(_MSC_VER)
    if (s_buf.st_mode & S_IFDIR) {
      return 1; // Directory
    } else if (s_buf.st_mode & S_IFREG) {
      return 2; // File
    } else {
      return 0; // Something else
    }
#else
    if (S_ISDIR(s_buf.st_mode) != 0) {
      return 1;
    } else if (S_ISREG(s_buf.st_mode) != 0) {
      return 2;
    } else {
      return 0;
    }
#endif
  } else { // sometimes files aren't properly determined so this is a fail safe
           // check
    fstream tmpObj(newPath.c_str(), ios::in);

    if (tmpObj.good()) {
      return 2;
      tmpObj.close();
    } else {
      //		printf("FILE: %s\n", newPath.c_str());
      //		cout << "errno: " << std::strerror(errno) << endl;
      return 0;
    }
  }
}

inline bool delete_path_tree(const char *path) {
  string tempStr(path);

  if ((strcmp(tempStr.c_str(), "") == 0) ||
      (strcmp(tempStr.c_str(), "NULL") == 0))
    return false;

  int pathAttrib = path_attribute(tempStr.c_str());

  if (pathAttrib == 1) // directory
  {
#if defined(_MSC_VER)
    string cmd = "rmdir /Q /S " + tempStr;
    CMD(cmd.c_str());
    return true;
#else
    string cmd = "rm -rf " + tempStr;
    CMD(cmd.c_str());
    return true;
#endif
  }

  if (pathAttrib == 2) // file
  {
#if defined(_MSC_VER)
    std::remove(tempStr.c_str());
    return true;
#else
    //		cout << "Deleting:\n" << tempStr << endl;
    string cmd = "rm -f " + tempStr;
    CMD(cmd.c_str());
    return true;
    return true;
#endif
  }

  return false;
}

inline bool do_mkdir(const char *path) {
#ifdef _MSC_VER
  string command = "mkdir " + (string)path;
  int result = path_attribute(path);
  if (result != 1)
    return (bool)system(command.c_str());
#else
  string command = "mkdir -p " + (string)path;
  int result = path_attribute(path);

  if (result != 1)
    return (bool)system(command.c_str());
#endif

  return false;
}

inline bool outFileExist(string fileName) {
  ifstream f(fileName.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }
}

inline void appendFileSuffix(char *oldFileWOExt, char *suffix,
                             string &newFileWOExt) {
  newFileWOExt = (string)oldFileWOExt + (string)suffix;
}

inline bool copyFile(char *source, char *destination) {
  string src(source);
  string dest(destination);

  if (strcmp(src.c_str(), dest.c_str()) == 0)
    return true;

  int rval;
#if defined(_MSC_VER) // implemented in case future implementation can be
                      // optimized
  delete_path_tree(dest.c_str());
  rval = CopyFileA(src.c_str(), dest.c_str(), 0);
#else
  //	cout << "Copying:\n" << src << endl << dest << endl;
  string command = "cp -f " + src + " " + dest;
  rval = CMD(command.c_str());
#endif

  return (bool)rval;
}

inline bool moveFile(
    char *source,
    char *destination) // NB: destination should be destination directory path
{
  string src = (string)source;
  string dest = (string)destination;

  string fullPath = src;

  std::size_t found = fullPath.find_last_of("//\\/");
  string fileName = fullPath.substr(found + 1);

  //	string newDest = dest + BKSLSH + fileName;
  string newDest = dest + fileName;

#if defined(_MSC_VER) // implemented in case future implementation can be
                      // optimized

  if (path_attribute(dest.c_str()) != 1)
    do_mkdir(dest.c_str());

  if (strcmp(newDest.c_str(), src.c_str()) == 0)
    return true;

  if (CopyFileA(src.c_str(), newDest.c_str(), false)) {
    delete_path_tree(src.c_str());
    return true;
  } else
    return false;
#else
  //	cout << "Moving:\n" << src << endl << newDest << endl;
  string command = "mv " + src + " " + newDest;
  if (CMD(command.c_str())) {
    return true;
  } else
    return false;
#endif
}

inline bool fileOperation(op command, string &source,
                          string destination = string(),
                          string suf = string()) {
  // enum{appendF, copyF, moveF, removeD, makeD, removeF} op;
  if (command == appendF) {
    appendFileSuffix((char *)destination.c_str(), (char *)suf.c_str(), source);
    return true;
  } else if (command == copyF) {
    return copyFile((char *)source.c_str(), (char *)destination.c_str());
  } else if (command == renameF) {
    if (copyFile((char *)source.c_str(), (char *)destination.c_str())) {
      return delete_path_tree(source.c_str());
    } else
      return false;
  } else if (command == moveF) {
    return moveFile((char *)source.c_str(), (char *)destination.c_str());
  } else if (command == removeF || command == removeD) {
    return delete_path_tree(source.c_str());
  } else if (command == makeD) {
    return do_mkdir(source.c_str());
  } else {
    cerr << "ERROR: Invalid option";
  }
  cout << "Common: FileIOMacros: inline bool fileOperation(op command: command"
       << (int)command << '\n';
  exit(0);
}
//=======================================================================================================
inline void suppressConsoleOutput() {
  std::cout.setstate(std::ios_base::failbit);
}
inline void restoreConsoleOutput() { std::cout.clear(); }

//=======================================================================================================
// Geared towards using a FILE* directly as an istream object
class cfilebuf : public std::streambuf {
  FILE *file;
  char c;
  int underflow() {
    int value = fgetc(this->file);
    if (value != EOF) {
      c = value;
      this->setg(&c, &c, &c + 1);
      return c;
    }
    return std::char_traits<char>::eof();
  }

public:
  cfilebuf(FILE *file) : file(file) {}
  // to own or not to own? ~cfilebuf() { fclose(this->file; }
};

class icfilestream : private virtual cfilebuf, public std::istream {
public:
  icfilestream(FILE *file)
      : cfilebuf(file), std::ios(static_cast<std::streambuf *>(this)),
        std::istream(static_cast<std::streambuf *>(this)) {}
};
//=======================================================================================================

#endif
