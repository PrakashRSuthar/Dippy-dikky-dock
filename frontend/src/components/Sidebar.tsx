import { useState, useCallback } from 'react';
import { NavLink } from 'react-router-dom';
import {
  Home, Clock, FolderOpen, Settings, Menu, X, Dna
} from 'lucide-react';

interface NavItem {
  to: string;
  label: string;
  icon: React.ReactNode;
}

const navItems: NavItem[] = [
  { to: '/', label: 'Dock', icon: <Home className="w-5 h-5" /> },
  { to: '/history', label: 'History', icon: <Clock className="w-5 h-5" /> },
  { to: '/runs', label: 'Runs', icon: <FolderOpen className="w-5 h-5" /> },
  { to: '/settings', label: 'Settings', icon: <Settings className="w-5 h-5" /> },
];

const Sidebar = () => {
  const [collapsed, setCollapsed] = useState(false);
  const [mobileOpen, setMobileOpen] = useState(false);

  const toggleMobile = useCallback(() => {
    setMobileOpen(prev => !prev);
  }, []);

  const closeMobile = useCallback(() => {
    setMobileOpen(false);
  }, []);

  const sidebarWidth = collapsed ? 'w-16' : 'w-60';

  return (
    <>
      {/* Mobile hamburger */}
      <button
        onClick={toggleMobile}
        className="fixed top-4 left-4 z-50 p-2 bg-gray-800 border border-gray-700 rounded-lg text-gray-300 hover:text-white hover:bg-gray-700 transition-colors md:hidden"
        aria-label="Toggle menu"
      >
        {mobileOpen ? <X className="w-5 h-5" /> : <Menu className="w-5 h-5" />}
      </button>

      {/* Mobile overlay */}
      {mobileOpen && (
        <div
          className="fixed inset-0 z-30 bg-black/60 md:hidden"
          onClick={closeMobile}
        />
      )}

      {/* Sidebar */}
      <aside
        className={`fixed top-0 left-0 z-40 h-screen bg-gray-900 border-r border-gray-700 flex flex-col transition-all duration-300 ease-in-out
          ${sidebarWidth}
          ${mobileOpen ? 'translate-x-0' : '-translate-x-full'}
          md:translate-x-0`}
      >
        {/* Logo */}
        <div className={`flex items-center gap-3 px-4 h-16 border-b border-gray-700 flex-shrink-0 ${collapsed ? 'justify-center' : ''}`}>
          <Dna className="w-6 h-6 text-blue-400 flex-shrink-0" />
          {!collapsed && (
            <span className="text-lg font-bold text-white whitespace-nowrap">DippyDock</span>
          )}
        </div>

        {/* Collapse toggle (desktop only) */}
        <button
          onClick={() => setCollapsed(prev => !prev)}
          className="hidden md:flex items-center justify-center h-8 border-b border-gray-700 text-gray-400 hover:text-white hover:bg-gray-800 transition-colors"
          aria-label={collapsed ? 'Expand sidebar' : 'Collapse sidebar'}
        >
          <Menu className={`w-4 h-4 transition-transform ${collapsed ? '' : 'rotate-90'}`} />
        </button>

        {/* Nav links */}
        <nav className="flex-1 py-4 space-y-1 overflow-y-auto">
          {navItems.map((item) => (
            <NavLink
              key={item.to}
              to={item.to}
              end={item.to === '/'}
              onClick={closeMobile}
              className={({ isActive }) =>
                `flex items-center gap-3 mx-2 px-3 py-2.5 rounded-lg transition-colors ${
                  collapsed ? 'justify-center' : ''
                } ${
                  isActive
                    ? 'bg-blue-600/20 text-blue-400 border border-blue-600/30'
                    : 'text-gray-400 hover:text-white hover:bg-gray-800 border border-transparent'
                }`
              }
            >
              <span className="flex-shrink-0">{item.icon}</span>
              {!collapsed && <span className="text-sm font-medium">{item.label}</span>}
            </NavLink>
          ))}
        </nav>

        {/* Footer */}
        <div className={`px-4 py-3 border-t border-gray-700 flex-shrink-0 ${collapsed ? 'hidden' : ''}`}>
          <p className="text-xs text-gray-500 text-center">v0.1.0</p>
        </div>
      </aside>

      {/* Spacer for main content */}
      <div className={`hidden md:block flex-shrink-0 transition-all duration-300 ${sidebarWidth}`} />
    </>
  );
};

export default Sidebar;
