MEXEXT=`mexext`

TARGET=mdlsAddNoise.$(MEXEXT) \
	mdlsMexTemplate.$(MEXEXT) \
	mdlsWriteLog.$(MEXEXT) 

MEX=mex

ifdef DEBUG
MEXOPT=-DDEBUG -g
else
MEXOPT=-O
endif

all: $(TARGET)


mdlsReconstruct.$(MEXEXT):mdlsReconstruct.c
	$(MEX) $(MEXOPT) $<

mdlsDeconstructFast.$(MEXEXT):mdlsDeconstructFast.cpp
	$(MEX) $(MEXOPT) $<

mdlsWriteLog.$(MEXEXT):mdlsWriteLog.c
	$(MEX) $(MEXOPT) $<

mdlsPGMRead.$(MEXEXT):mdlsPGMRead.c
	$(MEX) $(MEXOPT) $<

mdlsPGMWrite.$(MEXEXT):mdlsPGMWrite.c
	$(MEX) $(MEXOPT) $<

mdlsDeconstruct.$(MEXEXT):mdlsDeconstruct.cpp
	$(MEX) $(MEXOPT) $<

mdlsAddNoise.$(MEXEXT):mdlsAddNoise.cpp
	$(MEX) $(MEXOPT) $<

mdlsCanonicalTransform.$(MEXEXT):mdlsCanonicalTransform.cpp
	$(MEX) $(MEXOPT) $<

mdlsInitCanonicalTransform.$(MEXEXT):mdlsInitCanonicalTransform.cpp
	$(MEX) $(MEXOPT) $<

mdlsMexTemplate.$(MEXEXT):mdlsMexTemplate.cpp
	$(MEX) $(MEXOPT) $<

#mdlsSparsify.$(MEXEXT) : mdlsSparsify.cpp
#	$(MEX) $(MEXOPT) $<

mdlsImgDif.$(MEXEXT) : mdlsImgDif.c
	$(MEX) $(MEXOPT) $<

clean:
	-rm *.$(MEXEXT)
	-rm *.o
