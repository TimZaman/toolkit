#ifndef KEEPALIVE_H
#define KEEPALIVE_H

#include <QDebug>

class KeepAlive
{
public:

	KeepAlive(){
		qDebug() << "KeepAlive::KeepAlive()";
		KeepAliveMM();
	}

	~KeepAlive(){
		qDebug() << "KeepAlive::~KeepAlive()";
		KeepAliveDestructorMM();
	}

#ifdef __APPLE__ 
	void KeepAliveMM();
	void KeepAliveDestructorMM();
#else
	//?
#endif 

};

#endif