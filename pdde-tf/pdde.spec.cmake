Summary: Numerical continuation software for delay-differential equations
Name: pdde-cont
Version: @PACKAGE_VERSION@
Release: 1%{?dist}
License: GPL v2
Group: Application/Scientific
Source: pdde-cont-@PACKAGE_VERSION@.tar.bz2
URL: N/A
Vendor: Unknown
Packager: Robert Szalai <s z a l a i _at_ m m . b m e . h u>
Buildroot: %{_tmppath}/%{name}-%{version}-buildroot

%description
Numerical continuation software for delay-differential equations

%prep
%setup

%build
cmake -G "Unix Makefiles" -D CMAKE_INSTALL_PREFIX=/usr/ -D CMAKE_BUILD_TYPE=RELWITHDEBINFO .

%install
cd $RPM_BUILD_DIR/%{name}-%{version}
make
make DESTDIR=$RPM_BUILD_ROOT install

%files
%{_bindir}/*
%{_includedir}/*
%{_datadir}/*

%clean
rm -rf $RPM_BUILD_ROOT
