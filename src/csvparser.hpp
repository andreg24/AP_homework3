#ifndef CSVPARSER_HPP_
#define CSVPARSER_HPP_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <string>
#include <variant>
#include <optional>
#include <map>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#pragma GCC diagnostic pop


using namespace std;
using namespace boost::accumulators;

class CSVParser {
public:
	//CONSTRUCTOR
	CSVParser(const string& input_file);

	//PARSER function
	void read();
	
	// OPERATOR ()
	//Access to a single element of the dataset through parentesis
	const variant<optional <string>, optional<double>> operator()(const int row, const int col) const;

	//BASIC STATISTIC
	//The following functions allow to do some basic statistical operations on specified column of the dataset
	double mean_col(const size_t col_idx) const;
	double var_col(const size_t col_idx) const;
	double median_col(const size_t col_idx) const;
	double std_dev(const size_t col_idx) const;
	double covar(const size_t col_idx1, const size_t col_idx2) const;
	double correlation_analysis(size_t col_idx1, size_t col_idx2) const;

	//FREQUENCY COUNT
	//The function counts how many times each element appears in the specified column.
  map<string, int> countFrequency(const size_t col_idx) const;

	//SUMMARY
	//Compute frequency count for every element in every column. 
	//Compute mean, median, standard deviation, correlation analisis for every column of double.
	//Results saved in the specified file.
	void summary(const string& filename) const;

	//CLASSIFICATON
	//Search in the specified column for the wanted element.
	//The function prints the corresponding rows in a .txt file
	void classification(const string wanted, const size_t col_idx, const string& filename) const;


	//ITERATOR CLASS
	class ColIterator{
	private:
	const CSVParser& csvparser;
	size_t col;

	public:
		//CONSTRUCTOR
		ColIterator(const CSVParser& csvp, size_t size ) : csvparser(csvp), col(size){};

		//OPERATOR++
		ColIterator& operator++(){
			++col;
			return *this;
		}

		//OPERATOR!=
		bool operator!= ( const ColIterator& other) const{
			return col != other.col;
		}

		//OPERATOR*
		const variant<vector<optional<string>>, vector<optional<double>>>& operator*() const{
			return csvparser.dataset[col];
		}
	};

	//ITERATOR TO FIRST COLUMN
	ColIterator begin() const{
		return ColIterator(*this, 0);
	}

	//ITERATOR TO LAST COLUMN
	ColIterator end() const{
		return ColIterator(*this, dataset.size());
	}


	

private:
	const string input_file;			//csv file
	vector<variant<vector<optional<string>>, vector<optional<double>>>> dataset;
	vector<string> header;		//name of the columns
	size_t size;					//number of column
};

#endif
