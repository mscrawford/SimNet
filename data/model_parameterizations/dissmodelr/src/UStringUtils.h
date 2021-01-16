//---------------------------------------------------------------------------
#pragma once

#ifndef UStringUtilsH
#define UStringUtilsH
//---------------------------------------------------------------------------
#endif

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int StrToInt(const string &str)
{
 istringstream is(str);
 int i;
 is >> i;
 return i;
}

double StrToFloat(const string &str)
{
 istringstream is(str);
 double f;
 is >> f;
 return f;
}

template <class out_type, class in_value>
out_type convert(const in_value & t)
{
 stringstream stream;
 stream << t; // insert value to stream
 out_type result; // store conversion?s result here
 stream >> result; // write value to result
 return result;
}

/*
You use convert() like this:

double d;
string salary;
string s="12.56";
d=convert <double> (s); //d equals 12.56
salary=convert <string> (9000.0);//salary equals ?9000?

convert and to_string are from: (accessed 25 august 2003)
http://www.zdnet.com.au/builder/program/java/story/0,2000034779,20272027,00.htm
*/

template <class T>
void to_string(string & result, const T & t)
{
 ostringstream oss; // create a stream
 oss << t; // insert value to stream
 result=oss.str(); // extract value and write it to result
}
