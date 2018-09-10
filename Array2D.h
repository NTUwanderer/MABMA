#ifndef _ARRAY2D_H_
#define _ARRAY2D_H_

#include <cstdio>

template <class T>
class Array2D {
public:
    Array2D(): _array(NULL), returner(this) {}

	Array2D(int _h, int _w)
	: w(_w), h(_h), returner(this)
	{
	    init(_h,_w);
	}
    ~Array2D()
    {
    	if(_array != NULL)
    		delete[] _array;
	}

    class Proxy;
	
	void init(int _h, int _w)
	{
		w = _w;
		h = _h;
		_array = new T[w*h];
	}

    int getWidth()
    {
    	return w;
    }

    int getHeight()
    {
    	return h;
    }
    
    Proxy& operator[](T index)
    {
		returner.UpdateProxy(index);
		return returner;    	
    }

    void clone(const Array2D& s)
    {
    	w = s.w;
    	h = s.h;
    	delete[] _array;

    	_array = new T[w*h];
    	for(int i = 0; i < w*h; i++)
    	{
    		_array[i] = s._array[i];
    	}
    }
	
private:    
    T* _array;
    int w;
    int h;
	
	Proxy returner;
};

template <class T>
class Array2D<T>::Proxy {
public:
    Proxy() : a(0), ref(NULL) { }
    Proxy(int index, Array2D* _ref) : a(index), ref(_ref) { }
    Proxy(Array2D* _ref) : a(0), ref(_ref) { }
	
	void UpdateProxy(int index)
	{
		a = index;
	}

    T& operator[] (int index2) {
        return ref->_array[ref->w*this->a + index2];
    }
private:
    int a;
    Array2D* ref;
};



#endif
