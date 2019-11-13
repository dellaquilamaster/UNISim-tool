CXXFLAGS = -g

all:
	$(MAKE) root
	$(MAKE) -C main  ;			$(MAKE) -C main install

root:
	$(MAKE) -C UNISFramework ;		$(MAKE) -C UNISFramework install
	$(MAKE) -C UNISRootEvent ;		$(MAKE) -C UNISRootEvent install
	$(MAKE) -C detectors/DetectionSetup ;	$(MAKE) -C detectors/DetectionSetup install
	$(MAKE) -C detectors/Strip ;		$(MAKE) -C detectors/Strip install
	$(MAKE) -C detectors/Lamp ;		$(MAKE) -C detectors/Lamp install
	$(MAKE) -C generator ;			$(MAKE) -C generator install
	$(MAKE) -C UNISShared ;			$(MAKE) -C UNISShared install
	$(MAKE) -C LISETools ;			$(MAKE) -C LISETools install

install:
	$(MAKE) -C UNISFramework install
	$(MAKE) -C UNISRootEvent install
	$(MAKE) -C detectors/DetectionSetup install
	$(MAKE) -C detectors/Strip install
	$(MAKE) -C detectors/Lamp install
	$(MAKE) -C generator install
	$(MAKE) -C UNISShared install
	$(MAKE) -C LISETools install

distclean:
	$(MAKE) -C UNISFramework distclean
	$(MAKE) -C UNISRootEvent distclean
	$(MAKE) -C detectors/DetectionSetup distclean
	$(MAKE) -C detectors/Strip distclean
	$(MAKE) -C detectors/Lamp distclean
	$(MAKE) -C generator distclean
	$(MAKE) -C UNISShared distclean
	$(MAKE) -C LISETools distclean
	
	$(MAKE) -C main distclean

clean:
	$(MAKE) -C UNISFramework clean
	$(MAKE) -C UNISRootEvent clean
	$(MAKE) -C detectors/DetectionSetup clean
	$(MAKE) -C detectors/Strip clean
	$(MAKE) -C detectors/Lamp clean
	$(MAKE) -C generator clean
	$(MAKE) -C UNISShared clean
	$(MAKE) -C LISETools clean

	$(MAKE) -C main clean
