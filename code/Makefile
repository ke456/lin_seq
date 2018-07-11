
DIRS = sage_output magma_output ZZ_extra vec_ZZ_p_extra lzz_p_extra vec_lzz_p_extra lzz_pX_middle_product lzz_pX_extra lzz_pX_CRT ZZ_pX_extra ZZX_extra ZZ_CRT ZZ_pX_CRT mat_lzz_pX_extra mat_ZZ_p_extra mat_ZZ_pX_extra bivariate bivariate_transpose lin_seq

all:	
	$(foreach dir, $(DIRS), cd $(dir)/src ; make ; cd ../.. ;)

clean:
	rm -f lib/*
	rm -f include/*
	$(foreach dir, $(DIRS), cd $(dir)/src ; make clean ; cd ../.. ;)
	$(foreach dir, $(DIRS), cd $(dir)/test ; make clean ; cd ../.. ;)
