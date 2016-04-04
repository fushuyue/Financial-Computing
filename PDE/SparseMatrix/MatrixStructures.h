#include <iostream>
#include <iomanip>

#define DEBUG_PRINT true

struct SparseMatrix
{
public:
	struct element
	{
	public:
		float _value;
		int _i;
		int _j;
	};
	
	SparseMatrix()
	{
		_matrix = NULL;
		_n = 0;
		_m = 0;
		_len_element_list = 0;
	};
	
	element * _matrix;
	int _n, _m;
	int _len_element_list;
	
	/* part of your assigment is to make this work*/
	float get_value(int i, int j)
	{
		for(int ii = 0; ii < _len_element_list; ii++)
			if(_matrix[ii]._i == i && _matrix[ii]._j==j )
				return _matrix[ii]._value;
		return 0;
	}
	
	/* it is not required, but you could build this functions to make other work easier.*/
	void set_value(int i, int j , int index, float value)
	{
		_matrix[index]._i=i;
		_matrix[index]._j=j;
		_matrix[index]._value=value;
	}
	
	SparseMatrix operator *(SparseMatrix keke)
	{

		// build a m*n matrix to store the value
		float **a=new float*[_m];
		for(int i = 0; i < _m; ++i)
			a[i]=new float [_n];
		
		// initialize the matrix
		for(int i = 0; i < _m; ++i)
			for(int j = 0; j < _n; ++j)
				a[i][j]=0;
		
		SparseMatrix result;
		result._m=_m;
		result._n=_n;
		//std::cout<<keke.get_value(0,0)<<std::endl;
		for(int ii = 0; ii < this->_len_element_list; ii++)
			for(int jj=0; jj< keke._len_element_list; jj++)
			{
				float temp;
				if(this->_matrix[ii]._j==keke._matrix[jj]._i)
				{
					temp=this->_matrix[ii]._value*keke._matrix[jj]._value;
					a[this->_matrix[ii]._i][keke._matrix[jj]._j]+=temp;
				}
			}
	
		 //initialize the matrix
		
		int count=0;
		for(int i = 0; i < _m; ++i)
			for(int j = 0; j < _n; ++j)
			{
				if(a[i][j]!=0)
				{
					count+=1;
				}
			}
		result._matrix = new SparseMatrix::element[count];
		count=0;
		for(int i = 0; i < _m; ++i)
			for(int j = 0; j < _n; ++j)
			{
				if(a[i][j]!=0)
				{
					result.set_value(i, j, count, a[i][j]);
					count+=1;
				}
			}
		
		result._len_element_list=count;

		return result;
	}
	
	/*checks to make sure there are no duplicate elemens
	 and all i,j coordinates are between 0 and _n, _m respecivle*/
	bool is_valid()
	{
		for(int ii = 0; ii < _len_element_list; ii++)
		{
			if(_matrix[ii]._i >= _n)
	  {
		  if(DEBUG_PRINT)
			  std::cerr<<"_i coordinate out of range"<<std::endl;
		  return false;
	  }
			if(_matrix[ii]._j >= _m)
	  {
		  if(DEBUG_PRINT)
			  std::cerr<<"_j coordinate out of range"<<std::endl;
		  return false;
	  }
			
			for(int jj = ii+1; jj < _len_element_list; jj++)
				if(_matrix[ii]._i == _matrix[jj]._i && _matrix[ii]._j == _matrix[jj]._j)
				{
					if(DEBUG_PRINT)
					{
						std::cerr<<"duplicate element entry"<<std::endl;
					}
					return false;
				}
		}
		return true;
	}
	
	void print()
	{
		/* this is a bad way to work with sparse matricies */
		/*feel free to build out this sparse matrix API to make sure the elements are sorted*/
		int i;
		int j;
		bool not_exist = true;
		for(int ii = 0; ii < _n*_m; ii++)
		{
			i = int(float(ii)/float(_n));
			j = ii%_n;
			not_exist = true;
			for(int jj = 0; jj < _len_element_list; jj++)
				if(_matrix[jj]._i == i && _matrix[jj]._j == j)
				{
					std::cout << std::fixed << std::setw( 11 ) << std::setprecision( 4 ) << _matrix[jj]._value << "   ";
					not_exist = false;
					break;
				}
			
			if(not_exist)
				std::cout << std::fixed << std::setw( 11 ) << std::setprecision( 4 ) << 0.0 << "   ";
			
			if(j == _n - 1)
				std::cout<<'\n';
		}
	}
};
