// Filter system types - fully reusable for any filtering scenario

export interface FilterOption {
  label: string;
  value: string;
}

export interface FilterConfig {
  key: string;
  label: string;
  options: FilterOption[];
}

export interface SortOption {
  label: string;
  value: string;
}

export type SortDirection = 'asc' | 'desc';

export interface FilterState {
  filters: Record<string, string[]>;
  sort: string;
  sortDirection: SortDirection;
  searchQuery: string;
  searchField: string;
}

// Utility type for emitting combined sort value
export const getSortValue = (field: string, direction: SortDirection): string => {
  return `${field}_${direction}`;
};
