CXXFLAGS = -g

all:
	$(MAKE) root
	$(MAKE) -C main  ;			$(MAKE) -C main install

root:
	$(MAKE) -C UNISFramework ;		$(MAKE) -C UNISFramework install
	$(MAKE) -C UNISRootEvent ;		$(MAKE) -C UNISRootEvent install
	$(MAKE) -C UNISDetectedRootEvent ;	$(MAKE) -C UNISDetectedRootEvent install
	$(MAKE) -C detectors/DetectionSetup ;	$(MAKE) -C detectors/DetectionSetup install
	$(MAKE) -C detectors/Strip ;		$(MAKE) -C detectors/Strip install
	$(MAKE) -C detectors/Lamp ;		$(MAKE) -C detectors/Lamp install
	$(MAKE) -C detectors/Fazia ;		$(MAKE) -C detectors/Fazia install
	$(MAKE) -C detectors/Silicon ;		$(MAKE) -C detectors/Silicon install
	$(MAKE) -C detectors/Oscar ;		$(MAKE) -C detectors/Oscar install
	$(MAKE) -C generator ;			$(MAKE) -C generator install
	$(MAKE) -C UNISShared ;			$(MAKE) -C UNISShared install
	$(MAKE) -C UNISTarget ;		$(MAKE) -C UNISTarget install

install:
	$(MAKE) -C UNISFramework install
	$(MAKE) -C UNISRootEvent install
	$(MAKE) -C UNISDetectedRootEvent install
	$(MAKE) -C detectors/DetectionSetup install
	$(MAKE) -C detectors/Strip install
	$(MAKE) -C detectors/Lamp install
	$(MAKE) -C detectors/Fazia install
	$(MAKE) -C detectors/Silicon install
	$(MAKE) -C detectors/Oscar install
	$(MAKE) -C generator install
	$(MAKE) -C UNISShared install
	$(MAKE) -C UNISTarget install

distclean:
	$(MAKE) -C UNISFramework distclean
	$(MAKE) -C UNISRootEvent distclean
	$(MAKE) -C UNISDetectedRootEvent distclean
	$(MAKE) -C detectors/DetectionSetup distclean
	$(MAKE) -C detectors/Strip distclean
	$(MAKE) -C detectors/Lamp distclean
	$(MAKE) -C detectors/Fazia distclean
	$(MAKE) -C detectors/Silicon distclean
	$(MAKE) -C detectors/Oscar distclean
	$(MAKE) -C generator distclean
	$(MAKE) -C UNISShared distclean
	$(MAKE) -C UNISTarget distclean
	
	$(MAKE) -C main distclean

clean:
	$(MAKE) -C UNISFramework clean
	$(MAKE) -C UNISDetectedRootEvent clean
	$(MAKE) -C detectors/DetectionSetup clean
	$(MAKE) -C detectors/Strip clean
	$(MAKE) -C detectors/Lamp clean
	$(MAKE) -C detectors/Fazia clean
	$(MAKE) -C detectors/Silicon clean
	$(MAKE) -C detectors/Oscar clean
	$(MAKE) -C generator clean
	$(MAKE) -C UNISShared clean
	$(MAKE) -C UNISTarget clean

	$(MAKE) -C main clean
