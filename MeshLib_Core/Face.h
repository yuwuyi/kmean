#ifndef _FACE_H_
#define _FACE_H_

class Halfedge;

class Face
{
public:
	Face(){ m_halfedge = NULL; m_propertyIndex=-1;}
	~Face(){;}

	//Pointers for Halfedge Data Structure
	Halfedge    *	& he() { return m_halfedge; }
	Halfedge    *he() const { return m_halfedge; }
	//optional
	int				& index() {return m_propertyIndex; }
	std::string		& PropertyStr() { return m_propertyStr;}

protected:
	//for Halfedge Data Structure
	Halfedge	*	m_halfedge;

	//optional	
	std::string		m_propertyStr;
	int				m_propertyIndex; // index to Property array
};


#endif