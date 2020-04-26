#pragma once
#pragma warning(disable : 4996)
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <ctime> 
using namespace std;
class File_helper
{
public:
	File_helper() {}
	virtual ~File_helper() {}
	template<typename T>
	void write_vector_to_file(const vector<T>& myVector, string& filename)
	{
		ofstream ofs(filename, ios::out | ofstream::binary);
		ostream_iterator<char> osi{ ofs };
		const char* beginByte = (char*)&myVector[0];

		const char* endByte = (char*)&myVector.back() + sizeof(T);
		copy(beginByte, endByte, osi);
	}

	template<typename T>
	void read_vector_from_file(vector<T>& newVector, string& filename)
	{
		newVector.clear();
		vector<char> buffer{};
		ifstream ifs(filename, ios::in | ifstream::binary);
		istreambuf_iterator<char> iter(ifs);
		istreambuf_iterator<char> end{};
		copy(iter, end, back_inserter(buffer));
		newVector.resize(buffer.size() / sizeof(T));
		memcpy(&newVector[0], &buffer[0], buffer.size());
	}
	string get_time_str()
	{
		time_t rawtime;
		struct tm* timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer, sizeof(buffer), "%d-%m-%Y-%H_%M_%S", timeinfo);
		std::string str(buffer);
		return str;
	}

};

