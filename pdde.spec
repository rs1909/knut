Summary: Numerical continuation software for delay-differential equations
Name: pdde-cont
Version: 0.9.14
Release: 1
License: Custom
Group: Applications/Scientific
Source: pdde-cont-0.9.14.tar.gz
URL: N/A
Distribution: Fedora 4
Vendor: Unknown
Packager: Robert Szalai <s z a l a i _at_ m m . b m e . h u>
Buildroot: %{_tmppath}/%{name}-%{version}-buildroot

%description
Numerical continuation software for delay-differential equations

%prep
%setup

%build
%configure 

%install
cd $RPM_BUILD_DIR/%{name}-%{version}
make
make DESTDIR=${RPM_BUILD_ROOT} install

%files
%{_bindir}/*
%{_includedir}/*

%clean
rm -rf $RPM_BUILD_ROOT
