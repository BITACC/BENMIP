#pragma once
class BaseAdvancedVector
{
public:
	BaseAdvancedVector() :val_(0){}
	BaseAdvancedVector(int val) :val_(val){}
	int get(){ return val_; }
	void set(int val){ val_ = val; }
	
	~BaseAdvancedVector();

private:
	int val_;
	

};
