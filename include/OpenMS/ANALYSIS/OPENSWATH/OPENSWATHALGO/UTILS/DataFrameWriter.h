/*
 * DataFrameWriter.h
 *
 *  Created on: Aug 9, 2012
 *      Author: witold
 */

#ifndef DATAFRAMEWRITER_H_
#define DATAFRAMEWRITER_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

namespace OpenSwath
{
	struct IDataFrameWriter
	{
			virtual ~IDataFrameWriter(){}
			virtual void colnames(const std::vector<std::string> & colnames)=0;
			virtual void store(const std::string & rowname,
					const std::vector<double> & values) = 0;
	};

	struct DataMatrix: IDataFrameWriter
	{
		private:
			std::vector<std::string> rownames_;
			std::vector<std::vector<double> > store_;
			std::vector<std::string> colnames_;
		public:
			DataMatrix() :
					colnames_(), store_()
			{
			}
			void store(const std::string & rowname,
					const std::vector<double> & values)
			{
				rownames_.push_back(rowname);
				store_.push_back(values);
			}

			void colnames(const std::vector<std::string> & colnames){
				colnames_=colnames;
			}
	};

	struct CSVWriter: IDataFrameWriter
	{
			std::ofstream file_stream_;
			std::string sep_;
			std::string eol_;
			CSVWriter(std::string filename) :
					sep_("\t"), eol_("\n")
			{
				file_stream_.open(filename.c_str());
			}

			void store(const std::string & rowname,
					const std::vector<double> & values)
			{
				file_stream_ << rowname;
				file_stream_ << sep_;
				std::size_t ncol = values.size();
				for (size_t i = 0; i < ncol; ++i) {
					file_stream_ << std::setprecision(5) << values[i];
					if (i < (ncol - 1))
						file_stream_ << sep_;
				}
				file_stream_ << eol_; //append line-end
			}

			virtual ~CSVWriter(){
				file_stream_.flush();
				file_stream_.close();
				std::cout << "have flushed and closed the file stream" << std::endl;
			}

			void colnames(const std::vector<std::string> & colnames){
								std::size_t ncol = colnames.size();
								for (size_t i = 0; i < ncol; ++i) {
									file_stream_ << colnames[i];
									if (i < (ncol - 1))
										file_stream_ << sep_;
								}
								file_stream_ << eol_; //append line-end
			}
	};
}

#endif /* DATAFRAMEWRITER_H_ */
