MEXEXT=`mexext`

TARGET=\
	mdlsCanonicalTransform.$(MEXEXT) \
	mdlsInitCanonicalTransform.$(MEXEXT) \
	mdlsGetPatches.$(MEXEXT) \
	mdlsReconstruct.$(MEXEXT) \
	mdlsReconstructWithMask.$(MEXEXT) \
	mdlsDeconstruct.$(MEXEXT) \
	mdlsPGMRead.$(MEXEXT) \
	mdlsImgDif.$(MEXEXT)  \
	mdlsDeconstructFast.$(MEXEXT) \
	mdlsPGMWrite.$(MEXEXT) 

MEX=mex

ifdef DEBUG
MEXOPT=-DDEBUG -g
else
MEXOPT=-O
endif

all: $(TARGET)


mdlsReconstruct.$(MEXEXT):mdlsReconstruct.c
	$(MEX) $(MEXOPT) $<

mdlsReconstructWithMask.$(MEXEXT):mdlsReconstructWithMask.c
	$(MEX) $(MEXOPT) $<

mdlsGetPatches.$(MEXEXT):mdlsGetPatches.cpp
	$(MEX) $(MEXOPT) $<

mdlsDeconstructFast.$(MEXEXT):mdlsDeconstructFast.cpp
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

mdlsImgDif.$(MEXEXT) : mdlsImgDif.c
	$(MEX) $(MEXOPT) $<

clean:
	-rm *.$(MEXEXT)
	-rm *.o
