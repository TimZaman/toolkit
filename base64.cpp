

#include "base64.h"
#include <iostream>

#include <stdint.h>
#include <stdlib.h>


static char encoding_table[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                                'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
                                'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
                                'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
                                'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                                'w', 'x', 'y', 'z', '0', '1', '2', '3',
                                '4', '5', '6', '7', '8', '9', '+', '/'};
static int mod_table[] = {0, 2, 1};



void b64::build_decoding_table() {

	this->decoding_table = (char *) malloc(256);

	for (int i = 0; i < 64; i++){
		this->decoding_table[(unsigned char) encoding_table[i]] = i;
	}
}

void b64::base64_cleanup() {
	free(this->decoding_table);
}

b64::b64(){ //c'tor
	build_decoding_table();
}

b64::~b64(){ //d'tor
	base64_cleanup();
}




std::string b64::base64_encode(std::string src_str){
	//Public interface  
	unsigned int dst_len = 0;
	char * dst_data = base64_encode((const unsigned char *)src_str.c_str(), src_str.size(), &dst_len);

	std::string ret(dst_data, dst_len);
	delete dst_data;

	return ret;
}

std::string b64::base64_decode(std::string src_str){
	//Public interface
	unsigned int dst_len = 0;
	unsigned char * dst_data = base64_decode(src_str.c_str(), src_str.size(), &dst_len);

	std::string ret((char*)dst_data, dst_len);
	delete dst_data;

	return ret;
}

std::string b64::base64_encode(unsigned char const* src_data, unsigned int src_len){
	//Public interface

	unsigned int dst_len = 0;
	char * dst_data = base64_encode(src_data, src_len, &dst_len);

	std::string ret(dst_data, dst_len);
	delete dst_data;

	return ret;
}

std::string b64::base64_decode(char const* src_data, unsigned int src_len){
	//Public interface

	unsigned int dst_len = 0;
	unsigned char * dst_data = base64_decode(src_data, src_len, &dst_len);

	std::string ret((char*)dst_data, dst_len);
	delete dst_data;

	return ret;
}




char *  b64::base64_encode(const unsigned char *data,
                    unsigned int input_length,
                    unsigned int *output_length) {

	*output_length = 4 * ((input_length + 2) / 3);

	char *encoded_data = (char *)malloc(*output_length);
	if (encoded_data == NULL) return NULL;

	for (int i = 0, j = 0; i < input_length;) {
		uint32_t octet_a = i < input_length ? (unsigned char)data[i++] : 0;
		uint32_t octet_b = i < input_length ? (unsigned char)data[i++] : 0;
		uint32_t octet_c = i < input_length ? (unsigned char)data[i++] : 0;

		uint32_t triple = (octet_a << 0x10) + (octet_b << 0x08) + octet_c;

		encoded_data[j++] = encoding_table[(triple >> 3 * 6) & 0x3F];
		encoded_data[j++] = encoding_table[(triple >> 2 * 6) & 0x3F];
		encoded_data[j++] = encoding_table[(triple >> 1 * 6) & 0x3F];
		encoded_data[j++] = encoding_table[(triple >> 0 * 6) & 0x3F];
	}

	for (int i = 0; i < mod_table[input_length % 3]; i++){
		encoded_data[*output_length - 1 - i] = '=';
	}

	return encoded_data;
}


unsigned char * b64::base64_decode(const char *data,
                             unsigned int input_length,
                             unsigned int *output_length) {

	if (this->decoding_table == NULL) build_decoding_table();

	if (input_length % 4 != 0) return NULL;

	*output_length = input_length / 4 * 3;
	if (data[input_length - 1] == '=') (*output_length)--;
	if (data[input_length - 2] == '=') (*output_length)--;

	unsigned char *decoded_data = (unsigned char *) malloc(*output_length);
	if (decoded_data == NULL) return NULL;

	for (int i = 0, j = 0; i < input_length;) {
		uint32_t sextet_a = data[i] == '=' ? 0 & i++ : this->decoding_table[data[i++]];
		uint32_t sextet_b = data[i] == '=' ? 0 & i++ : this->decoding_table[data[i++]];
		uint32_t sextet_c = data[i] == '=' ? 0 & i++ : this->decoding_table[data[i++]];
		uint32_t sextet_d = data[i] == '=' ? 0 & i++ : this->decoding_table[data[i++]];

		uint32_t triple = (sextet_a << 3 * 6)
		+ (sextet_b << 2 * 6)
		+ (sextet_c << 1 * 6)
		+ (sextet_d << 0 * 6);

		if (j < *output_length) decoded_data[j++] = (triple >> 2 * 8) & 0xFF;
		if (j < *output_length) decoded_data[j++] = (triple >> 1 * 8) & 0xFF;
		if (j < *output_length) decoded_data[j++] = (triple >> 0 * 8) & 0xFF;
	}

	return decoded_data;
}





