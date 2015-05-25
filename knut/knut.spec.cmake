Summary: Numerical continuation software for delay-differential equations
Name: Knut
Version: @PACKAGE_VERSION@
Release: 1%{?dist}
License: GPL v2
Group: Application/Scientific
Source: Knut-@PACKAGE_VERSION@.tar.bz2
URL: https://gitorious.org/knut
Vendor: Unknown
Packager: Robert Szalai <robicjedi _at_ gmail.com>
Buildroot: %{_tmppath}/%{name}-%{version}-buildroot
Buildrequires: cmake >= 2.8.11
Buildrequires: qt5-qtbase-devel >= 5.1
Buildrequires: qt5-qtbase-gui >= 5.1
Buildrequires: qt5-qtsvg-devel >= 5.1
Buildrequires: openblas-devel
Buildrequires: suitesparse-devel
Buildrequires: mxml-devel >= 2.8
Buildrequires: gcc-gfortran >= 4.7.2
Buildrequires: libgfortran-static >= 4.7.2
Buildrequires: gcc-c++ >= 4.7.2
Buildrequires: libstdc++-devel >= 4.7.2
Buildrequires: libgfortran >= 4.7.2

%description
Numerical continuation software for delay-differential equations

%prep
%setup

%build
%cmake
make %{?_smp_mflags}

%install
cd $RPM_BUILD_DIR/%{name}-%{version}
make DESTDIR=$RPM_BUILD_ROOT install

%files
%{_bindir}/*
%{_includedir}/*
%{_datadir}/*

%clean
rm -rf $RPM_BUILD_ROOT

%changelog
* Tue Dec 31 2013 Robert Szalai <robicjedi@gmail.com> 9-1
- The first rpm release
