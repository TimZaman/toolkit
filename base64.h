#ifndef BASE64_TIM_H
#define BASE64_TIM_H

#include <string>



class b64 {
	 
	public:
		b64();
		~b64();
		std::string base64_encode(unsigned char const* , unsigned int len);
		std::string base64_encode(std::string);
		std::string base64_decode(char const* , unsigned int len);
		std::string base64_decode(std::string);
	private:
		char *decoding_table = NULL;

		void base64_cleanup();
		void build_decoding_table();

		char * base64_encode(const unsigned char *data, unsigned int input_length, unsigned int *output_length);
		unsigned char * base64_decode(const char *data, unsigned int input_length, unsigned int *output_length);

};

#endif