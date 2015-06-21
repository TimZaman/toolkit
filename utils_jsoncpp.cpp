//utils_jsoncpp.cpp

#include "utils_jsoncpp.h"





template <typename Iterable>
Json::Value util::iterable2json(Iterable const& cont) {
    Json::Value v;
    for (auto&& element: cont) {
        v.append(element);
     }
     return v;
}

void util::PrintJSONValue( Json::Value val ){
	if( val.isString() ) {
		printf( "string(%s)", val.asString().c_str() ); 
	} else if( val.isBool() ) {
		printf( "bool(%d)", val.asBool() ); 
	} else if( val.isInt() ) {
		printf( "int(%d)", val.asInt() ); 
	} else if( val.isUInt() ) {
		printf( "uint(%u)", val.asUInt() ); 
	} else if( val.isDouble() ) {
		printf( "double(%f)", val.asDouble() ); 
	} else {
		printf( "unknown type=[%d]", val.type() ); 
	}
}

bool util::PrintJSONTree( Json::Value &root, unsigned short depth /* = 0 */) {
	depth += 1;
	printf( " {type=[%d], size=%d}", root.type(), root.size() ); 

	if( root.size() > 0 ) {
		printf("\n");
		for( Json::ValueIterator itr = root.begin() ; itr != root.end() ; itr++ ) {
			// Print depth. 
			for( int tab = 0 ; tab < depth; tab++) {
			   printf("-"); 
			}
			printf(" subvalue(");
			PrintJSONValue(itr.key());
			printf(") -");
			PrintJSONTree( *itr, depth); 
		}
		return true;
	} else {
		printf(" ");
		PrintJSONValue(root);
		printf( "\n" ); 
	}
	return true;
}

